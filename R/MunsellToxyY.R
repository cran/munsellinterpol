

#   MunsellSpec     a character vector of Munsell names
#                   or a numeric matrix with 3 columns, with HVC in each row
#                   or a single plain numeric HVC 3-vector
#   xyC             xy 2-vector for illuminant C
#                   can also be a string: 'NBS', 'JOSA', 'NTSC', or 'ASTM'
#   hcinterp        'bicubic' or 'bilinear'
#   vinterp         'cubic' or 'linear'
#   space           'xy' or 'ab'
#   YfromV          'ASTM' or 'NBS' or 'MUNSELL'. 'MGO' is not allowed here

#   return value    a data.frame

MunsellToxyY  <-  function( MunsellSpec,
                            xyC='NBS',
                            hcinterp='bicubic',
                            vinterp='cubic',
                            YfromV='ASTM',
                            warn=TRUE )
    {
    time_start  = gettime()

    if( is.character(MunsellSpec) )
        {
        HVC         = HVCfromMunsellName(MunsellSpec)
        SAMPLE_NAME = MunsellSpec
        }
    else
        {
        HVC = prepareNx3( MunsellSpec )
        if( is.null(HVC) )  return(NULL)

        SAMPLE_NAME = MunsellNameFromHVC( HVC, digits=3 )
        #rownames(HVC)   = MunsellNameFromHVC( HVC, digits=3 )

        if( is.null(colnames(HVC)) )
            colnames(HVC)   = c("H","V","C")
        }

    #   process xyC
    xyC_org = xyC
    xyC     = process_xyC( xyC )
    if( any(is.na(xyC)) )
        {
        event_level( ERROR, "Argument xyC is invalid.", class="invalid_argument", extra=list(xyC=xyC_org) )
        return(NULL)
        }

    # process hcinterp
    hcinterp_org    = hcinterp
    hcinterp        = process_hcinterp(hcinterp)
    if( is.na(hcinterp) )
        {
        event_level( ERROR, "Argument hcinterp='%s' is invalid.", hcinterp_org,
                                class="invalid_argument", extra=list(hcinterp=hcinterp_org) )
        return(NULL)
        }

    #   assign function pointer
    if( hcinterp == 'bicubic' )
        hcinterpfun = bicubicCardinal   # mybicubic
    else if( hcinterp == 'bilinear' )
        hcinterpfun = bilinearApprox    # mybilinear


    # process vinterp
    vinterp_org = vinterp
    vinterp     = process_vinterp(vinterp)
    if( is.na(vinterp) )
        {
        event_level( ERROR, "Argument vinterp='%s' is invalid.", vinterp_org,
                            class="invalid_argument", extra=list(vinterp=vinterp_org) )
        return(NULL)
        }

    # process YfromV
    full    = c( 'ASTM', 'NBS', 'MUNSELL', 'PRIEST' )
    idx     = pmatch( toupper(YfromV), full )
    if( is.na(idx) )
        {
        event_level( ERROR, "Argument YfromV='%s' is invalid. ('MgO' is not allowed in this function).", YfromV,
                            class="invalid_argument", extra=list(YfromV=YfromV) )
        return(NULL)
        }

    #   ready to process data

    #   test any element of the lookup tables to see whether illuminant C has changed
    same    =  xyC[1]==p.LookupList[[1]]$x[1,1]  &&  xyC[2]==p.LookupList[[1]]$y[1,1]
    if( is.na(same)  ||  ! same )
        {
        base::unlockBinding( "p.LookupList", asNamespace('munsellinterpol') )

        # fill all the lookup tables with x and y for illuminant C
        log_level( TRACE, "filling p.LookupList with xy=%g,%g for illuminant C.", xyC[1], xyC[2] )
        for( k in 1:length(p.LookupList) )
            {
            p.LookupList[[k]]$x[ ,1] <<- xyC[1]
            p.LookupList[[k]]$y[ ,1] <<- xyC[2]
            }
        }

    n   = nrow(HVC)

    xyY = matrix( NA_real_, n, 3 )
    colnames( xyY ) = c( 'x', 'y', 'Y' )

    #   do all Y coordinates in only 1 call
    xyY[ ,3]    = YfromV( HVC[ ,2], which=YfromV  )  #*  (whiteC[3]/100)


    V.vector    = attr( p.LookupList, "V.vector" )
    H.vector    = attr( p.LookupList, "H.vector" )

    for( i in 1:n )
        {
        chroma  = HVC[i,3]
        if( ! is.finite(chroma) || chroma<0 )
            {
            #   log_level( ERROR, "Failed to map sample %d, of %d. chroma=%g.", i, n, chroma )
            next
            }

        if( chroma == 0 )
            {
            #   neutrals are a special case
            xyY[i,1:2] = xyC[1:2]
            next
            }

        hue     = HVC[i,1]
        value   = HVC[i,2]

        if( value == 0 )
            {
            #   we must have 0<chroma, and this is physically unrealizable
            next
            }

        iV = match( value, V.vector )
        if( ! is.na(iV) )
            {
            C.vector    = attr( p.LookupList[[iV]], "C.vector" )

            #   special case, the input point is exactly on one of the constant Value planes, so only 1 HC interpolation required.
            #   But this will probably happen a lot
            xyY[i,1]    = hcinterpfun( H.vector, C.vector, p.LookupList[[iV]]$x, hue, chroma )$z
            xyY[i,2]    = hcinterpfun( H.vector, C.vector, p.LookupList[[iV]]$y, hue, chroma )$z

            #if( any( is.na(xyY[i, ]) ) )
            #    log_level( ERROR, "Failed to map sample %d, of %d. Exactly on V plane.", i, n )

            next
            }

        vinterp.i   = vinterp

        if( p.vinterpOverride  &&  value < 2 )
            vinterp.i   = 'linear'

        #   value is *between* V-planes, so a lot more work is needed
        iV  = findInterval( value, V.vector, all.inside=TRUE )
        if( vinterp.i == 'linear' )
            vlevel = iV : (iV+1)
        else if( vinterp.i == 'cubic' )
            vlevel = max(iV-1,1) : min(iV+2, length(V.vector))

        xymat   = matrix( NA, length(vlevel), 2 )
        for( k in 1:length(vlevel) )
            {
            iV  = vlevel[k]

            C.vector    = attr( p.LookupList[[iV]], "C.vector" )

            xymat[k,1]    = hcinterpfun( H.vector, C.vector, p.LookupList[[iV]]$x, hue, chroma )$z
            xymat[k,2]    = hcinterpfun( H.vector, C.vector, p.LookupList[[iV]]$y, hue, chroma )$z
            }

        vvec    = V.vector[vlevel]

        if( any( is.na(xymat) ) )
            {
            #log_level( ERROR, "Failed to map sample %d, of %d. Between V planes.", i, n )
            #rownames( xymat ) = as.character(vvec)          #; colnames(xymat) = c('x','y') ;       print( xymat )
            #print( xymat )
            next
            }

        #   now interpolate along each column of xymat
        if( vinterp.i == 'linear' )
            {
            xyY[i,1]    = stats::approx( vvec, xymat[ ,1], xout=value )$y
            xyY[i,2]    = stats::approx( vvec, xymat[ ,2], xout=value )$y
            }
        else if( vinterp.i == 'cubic' )
            {
            #xyY[i,1]    = spline( vvec, xymat[ ,1], xout=value, method='natural' )$y
            #xyY[i,2]    = spline( vvec, xymat[ ,2], xout=value, method='natural' )$y

            xyY[i,1:2]    = splineCardinal( vvec, xymat, xout=value )
            }
        }


    out = data.frame( row.names=1:n, stringsAsFactors=FALSE )

    out$SAMPLE_NAME = SAMPLE_NAME
    out$HVC         = HVC

    rownames(xyY)   = rownames(HVC)
    if( is.null(rownames(xyY) ) )   rownames(xyY)   = SAMPLE_NAME
    out$xyY         = xyY

    if( FALSE  &&   100 <= n )
        {
        time_elapsed = gettime() - time_start
        log_level( INFO, "Processed %d samples in %g seconds.  %g/sample.",
                            n, time_elapsed, time_elapsed/n )
        }


    if( warn )
        {
        #   check for mapping failure
        bad = is.na(xyY[ ,1])  |  is.na(xyY[ ,2])

        if( any(bad) )
            {
            event_level( WARN, "%d samples, out of %d, could not be mapped; xy set to NA.", sum(bad), length(bad),
                                class="incomplete", extra=list( invalid=HVC[bad, ,drop=FALSE], indexes=which(bad) )  )
            }

        #   check for xy outside the triangle
        outside = xyY[ ,1]<0  |  xyY[ ,2]<0  |  1<xyY[ ,1]+xyY[ ,2]
        outside[ is.na(outside) ] = FALSE

        if( any(outside) )
            {
            res = event_level( WARN, "%d samples, out of %d, were mapped to xy outside the triangle.", sum(outside), n,
                                    class="incomplete", extra=list( invalid=HVC[outside, ,drop=FALSE], indexes=which(outside) ) )

            #if( ! is.null(res)  &&  sum(outside) <= 10 )
            #    {
            #    #  since res is not NULL, a logging event actually occurred, and the number of problems is small
            #    print( out[outside, ] )
            #    }
            }
        }

    #   print( str(out) )



    return(out)
    }

#   xyC     single string or numeric 2-vector
#
#   returns valid numeric 2-vector, or NA in case that xyC is invalid

process_xyC <- function( xyC )
    {
    if( is.character(xyC) )
        {
        idx = pmatch( toupper(xyC), rownames(p.xyC) )

        if( is.na(idx) )
            return(NA)

        xyC   = p.xyC[ idx, ]
        }

    ok  = is.numeric(xyC)  &&  length(xyC) == 2
    if( ! ok )
        return(NA)


    #   if( length(whiteC) == 2 )    whiteC = c( whiteC, 100 )

    return( xyC )
    }



process_hcinterp <- function( hcinterp )
    {
    full    = c('bicubic','bilinear')

    idx = pmatch( tolower(hcinterp), full )
    if( is.na(idx) )
        return(NA)

    return( full[idx] )
    }

process_vinterp <- function( vinterp )
    {
    full   = c('cubic','linear')
    idx = pmatch( tolower(vinterp), full )
    if( is.na(idx) )
        return(NA)

    return( full[idx] )
    }

