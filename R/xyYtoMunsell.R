#
#   xyY             a numeric Nx3 matrix, or a vector that can be converted to one,
#                   with xyY in each row
#   xyC             xy  2-vector for illuminant C
#                   can also be a string: 'NBS', 'JOSA', 'NTSC', or 'ASTM'
#   hcinterp        'bicubic' or 'bilinear'
#   vinterp         'cubic' or 'linear'
#   VfromY          'ASTM' or 'NBS' or 'MUNSELL'.  'MGO' is not allowed here
#   warn            log a warning if there were any inversion failures
#   perf            add performance data to returned output
#
#   return value    a data.frame with columns
#                   xyY             the input
#                   HVC             a matrix with the same size as xyY, with (Hue,Value,Chroma) in the rows
#                   SAMPLE_NAME     Munsell notation for HVC, a string

xyYtoMunsell  <-  function( xyY,
                            xyC='NBS',
                            hcinterp='bicubic',
                            vinterp='cubic',
                            VfromY='ASTM',
                            warn=TRUE,
                            perf=FALSE )
    {
    p   = 'rootSolve'
    if( ! requireNamespace( p, quietly=TRUE ) )
        {
        log.string( ERROR, "required package '%s' could not be loaded.", p )
        return(NULL)
        }

    time_start  = gettime()

    xyY = prepareNx3(xyY)
    if( is.null(xyY) )  return(NULL)
    
    if( is.null(colnames(xyY)) )
        colnames(xyY)   = c( 'x', 'y', 'Y' )

    #   process xyC
    xyC  = process_xyC( xyC )
    if( any(is.na(xyC)) )
        {
        log.string( ERROR, "xyC is invalid." )
        return(NULL)
        }

    # process hcinterp
    hcinterp    = process_hcinterp(hcinterp)
    if( is.na(hcinterp) )
        {
        log.string( ERROR, "hcinterp='%s' is invalid.", hcinterp )
        return(NULL)
        }

    #   assign function pointer
    if( hcinterp == 'bicubic' )
        hcinterpfun = bicubicCardinal  # mybicubic
    else if( hcinterp == 'bilinear' )
        hcinterpfun = bilinearApprox    # mybilinear



    # process vinterp
    vinterp    = process_vinterp(vinterp)
    if( is.na(vinterp) )
        {
        log.string( ERROR, "vinterp='%s' is invalid.", as.character(vinterp) )
        return(NULL)
        }


    #   process VfromY
    full    = c( 'ASTM', 'NBS', 'MUNSELL', 'PRIEST' )
    idx     = pmatch( toupper(VfromY), full )
    if( is.na(idx) )
        {
        log.string( ERROR, "VfromY='%s' is invalid. ('MgO' is not allowed in this function).", VfromY )
        return(NULL)
        }

    # allocate arrays for mapping H,C -> x and H,C -> y
    # Value 7 has the full extent, so just use that as a template
    xlookup = array( NA_real_, dim=dim(p.LookupList[[7]]$x) )
    ylookup = array( NA_real_, dim=dim(p.LookupList[[7]]$y) )

    H.vector    = attr( p.LookupList, "H.vector" )   # ; cat( 'H.vector=',H.vector,'\n')
    V.vector    = attr( p.LookupList, "V.vector" )
    C.vector    = attr( p.LookupList[[7]], "C.vector" )

    #   fill neutrals at Chroma=0, these are the same for all hues and values
    xlookup[ ,1] = xyC[1]
    ylookup[ ,1] = xyC[2]



    #   define the "forward" function, which we will try to find the root of
    #   AB          input 2-vector.  It changes on each iteration.
    #   value       the Munsell Value.  It remains constant throughout the iteration
    #   vlevel      vector of value indices (not the values) where HC is interpolated.
    #               This is 1-based  It depends only on value. It remains constant throughout the iteration
    #   xy_target   we are trying to find HC that maps to xy_target.  It remains constant throughout the iteration
    forwardfun <- function( AB,    value, vlevel, xy_target )
        {
        evaluation <<- evaluation + 1
        #cat( "######  evaluation = ", evaluation, " ###########\n", file=stderr() )
        # cat( 'HC=', HC, "  value=", value, '\n' )

        HC  = unlist( HCfromAB( AB[1], AB[2] ) )

        #HC[1] = HC[1] %% 100
        #
        #if( HC[2] < 0 )
        #    log.string( WARN, "Chroma = %g < 0.  hue=%g value=%g.  evaluation=%d.",
        #                HC[2], HC[1], value, evaluation )


        iH  = findInterval( HC[1], H.vector, all.inside=TRUE )  #; cat( 'iH=',iH,'\n')
        iC  = findInterval( HC[2], C.vector, all.inside=TRUE )

        if( hcinterp == 'bicubic' )
            {
            hseq    = max(iH-1,1) : min(iH+2,length(H.vector))      # 4 points for bicubic
            cseq    = max(iC-1,1) : min(iC+2,length(C.vector))      # 4 points for bicubic
            }
        else
            {
            hseq    = iH : (iH+1)      # 2 points for bilinear
            cseq    = iC : (iC+1)      # 2 points for bilinear
            }

        #   iterate over 4x4 or 2x2 grid, and compute missing values
        #print(hseq)
        #print(cseq)

        gridpoints = 0

        #   extract just the values that are relevant to the iteration, usually either 2 or 4 of them
        #   V.vector will not be used again
        Vsub    = V.vector[ vlevel ]

        for( ih in hseq )
            {
            for( ic in cseq )
                {
                if( ! is.na( xlookup[ih,ic] ) ) next

                #   xlookup[ih,ic] has been computed yet, so compute it now
                gridpoints = gridpoints + 1

                xvec    = rep( NA_real_, length(Vsub) )  # length(V.vector)
                yvec    = xvec

                for( k in 1:length(vlevel) )      # v in vlevel )
                    {
                    iV   = vlevel[k]
                    icc  = min( ic, ncol(p.LookupList[[iV]]$x) )
                    xvec[k] = p.LookupList[[iV]]$x[ih,icc]
                    yvec[k] = p.LookupList[[iV]]$y[ih,icc]
                    }

                fmask   = is.finite(xvec)
                count   = sum( fmask )    #; print(count)

                if( vinterp == 'linear'   ||  (p.vinterpOverride && value<2)  )
                    {
                    # print( xvec ) ; print( V.vector ) ; print( value )
                    if( count < 1 )
                        {
                        #   need 1 finite value
                        if( evaluation == 1 )
                            mess = sprintf( "forwardfun(). The initial guess for root is outside the range of the HVC -> xyY lookup table.  AB=%g,%g   HC=%g,%g",
                                                 AB[1], AB[2], HC[1], HC[2] )
                        else
                            mess = sprintf( "forwardfun(). The root-finder has wandered outside the range of the HVC -> xyY lookup table. evaluation=%d. AB=%g,%g   HC=%g,%g",
                                                evaluation, AB[1], AB[2], HC[1], HC[2] )
                        stop( mess, call.=FALSE )
                        }

                    if( count == 2 )  #  all(fmask)
                        {
                        #   interpolate between those 2
                        xlookup[ih,ic]  <<- approx( Vsub, xvec, xout=value )$y
                        ylookup[ih,ic]  <<- approx( Vsub, yvec, xout=value )$y
                        }
                    else
                        {
                        #   just use the one
                        #   this happens for 1 out of 994 optimals
                        idx = which( fmask )
                        xlookup[ih,ic]  <<- xvec[ idx ]
                        ylookup[ih,ic]  <<- yvec[ idx ]
                        }
                    }
                else if( vinterp == 'cubic' )
                    {
                    if( count < 2 )
                        {
                        #   need 2 finite values
                        if( evaluation == 1 )
                            mess = sprintf( "forwardfun(). The initial guess for root is outside the range of the HVC -> xyY lookup table.   AB=%g,%g   HC=%g,%g",
                                                AB[1], AB[2], HC[1], HC[2]  )
                        else
                            mess = sprintf( "forwardfun(). The root-finder has wandered outside the range of the HVC -> xyY lookup table. evaluation=%d.  AB=%g,%g  HC=%g,%g",
                                                evaluation, AB[1], AB[2], HC[1], HC[2] )
                        stop( mess, call.=FALSE )
                        }

                    xyvec =  cbind(xvec,yvec)
                    rownames(xyvec) = as.character(Vsub)
                    #print( xyvec )

                    if( all(fmask) )
                        #   the usual case
                        xy  = splineCardinal( Vsub, xyvec, xout=value )
                    else
                        {
                        #   despite the dilation, count==3 still happens for 11 out of 994 optimals
                        #   TODO:  another workaround might be shrinking the polym() model, as already done for value=0.2
                        # log.string( TRACE, "cubic  count=%d   evaluation=%d   value=%g", count, evaluation, value )
                        xy  = splineCardinal( Vsub[fmask], xyvec[fmask, ], xout=value )
                        }

                    #xlookup[ih,ic]  <<- spline( Vsub, xvec, xout=value, method='natural' )$y
                    #ylookup[ih,ic]  <<- spline( Vsub, yvec, xout=value, method='natural' )$y

                    xlookup[ih,ic]  <<- xy[1]
                    ylookup[ih,ic]  <<- xy[2]
                    }
                }
            }


        #   the grids xlookup and ylookup are ready, now interpolate in this V-plane
        xy      = numeric(2)
        xy[1]   = hcinterpfun( H.vector, C.vector, xlookup, HC[1], HC[2] )$z
        xy[2]   = hcinterpfun( H.vector, C.vector, ylookup, HC[1], HC[2] )$z


        #cat( 'vlevel=', vlevel, '\n' )
        #cat( 'xy_target=', xy_target, '\n' )
        #cat( 'gridpoints=', gridpoints, '\n' )
        #cat( 'xy=', xy, '\n' )

        delta   =  xy - xy_target
        #cat( 'delta=', sum(abs(delta)), '\n' )

        return(delta)
        }



    HVC = array( NA_real_, dim=dim(xyY) )
    colnames(HVC)   = c( 'H', 'V', 'C' )

    HVC[ ,2]    = VfromY( xyY[ ,3], which=VfromY )

    n   = nrow(xyY) #; print(n)

    evaluations.vec     = rep( NA_integer_, n )
    iterations.vec      = rep( NA_integer_, n )
    estim.precis.vec    = rep( NA_real_, n )
    time.elapsed.vec    = rep( NA_real_, n )

    for( i in 1:n )
        {
        time_begin  = gettime()  # as.double( Sys.time() )

        xy_target   = xyY[i,1:2]
        value       = HVC[i,2]

        if( is.na( value ) )    next

        if( 0 < value  &&  any( is.na(xy_target) ) )    next    # invalid, so ignore this one

        if( value==0  ||  all( xy_target == xyC ) )
            {
            #   neutrals are easy.  Note that value==0 must be pure black !
            #   Set both H and C to 0
            HVC[ i, c(1,3) ] = 0
            time.elapsed.vec[i] = (gettime()  - time_begin)
            next
            }

        delta   = abs(V.vector - value)     #; print( delta )
        iV  = which.min( delta )
        if( delta[iV] < 5.e-9 )
            {
            #   this is so close to a V-plane, that this is probably where xyY[i, ] came from
            value       = V.vector[iV]
            HVC[i,2]    = value
            }

        if( value < 0  ||  10 < value )
            #   this is out of range
            next

        iV  = match(value,V.vector)

        if( is.na(iV) )
            {
            #   value is *between* 2 V-planes. This is the usual case
            #   erase the lookup tables, but only where Chroma > 0
            xlookup[ ,-1] = NA_real_    #; print(xlookup)
            ylookup[ ,-1] = NA_real_

            iV  = findInterval( value, V.vector, all.inside=TRUE )
            if( vinterp == 'linear'  || (p.vinterpOverride && value<2) )
                vlevel = iV : (iV+1)
            else if( vinterp == 'cubic' )
                vlevel = max(iV-1,1) : min(iV+2, length(V.vector))
            }
        else
            {
            #   value is *exactly* on a V-plane
            #   not likely but since all the action will be on this plane, we can optimize it
            #   cat( "i=", i, "  Exact value=", value, '\n', file=stderr() )
            nC  = dim(p.LookupList[[iV]]$x)[2]
            xlookup[ ,1:nC] = p.LookupList[[iV]]$x
            ylookup[ ,1:nC] = p.LookupList[[iV]]$y
            vlevel  = iV    # vlevel will never be accessed, but this is what it would be if it were
            }

        #   make an initial guess for HC
        evaluation = 0

        #HCstart  = approxHCfromxyY( xyY[i, ] , c(xyC,100) )
        #ABstart = unlist( ABfromHC( HCstart[1], HCstart[2] ) )

        ABstart = approxABfromxyY( xyY[i, ], c(xyC,100), V.vector )
        HCstart = unlist( HCfromAB( ABstart[1], ABstart[2] ) )

        res = try( rootSolve::multiroot( forwardfun, ABstart, rtol=1.e-8, atol=1.e-6, ctol=0, verbose=FALSE,
                                            value=value, vlevel=vlevel, xy_target=xy_target ),  silent=F )

        if( inherits(res,"try-error" ) )   # class(res) == "try-error" )
           {
           #    failed to find the root
           # print(res)
           log.string( DEBUG, "rootSolve::multiroot() failed.  ABstart=%g,%g.   HCstart=%g,%g.  value=%g",
                                         ABstart[1], ABstart[2], HCstart[1], HCstart[2], value )
           next
           }

        time.elapsed.vec[i] = (gettime() - time_begin)    # ;as.double( Sys.time() ) - time_begin

        #   convert from AB back to HC
        HVC[ i, c(1,3) ] = unlist( HCfromAB( res$root[1], res$root[2] ) )

        #   cat( "Success !  iterations=", res$iter, '\n', file=stderr() )

        iterations.vec[i]   = res$iter
        evaluations.vec[i]  = evaluation
        estim.precis.vec[i] = res$estim.precis
        }


    #   hue wrap-around to interval (0,100]
    hue             = HVC[ ,1] %% 100
    hue[ hue==0 ]   = 100
    HVC[ , 1]       = hue

    SAMPLE_NAME = MunsellNameFromHVC(HVC)

    #   copy rownames from xyY to HVC, unless NULL then use SAMPLE_NAME
    rownames(HVC)   = rownames(xyY)
    if( is.null(rownames(HVC)) )    rownames(HVC) = SAMPLE_NAME

    out = data.frame( row.names=1:n, stringsAsFactors=FALSE )

    out$xyY         = xyY
    out$HVC         = HVC
    out$SAMPLE_NAME = SAMPLE_NAME

    if( perf )
        {
        out$time.elapsed    = time.elapsed.vec
        out$iterations      = iterations.vec
        out$evaluations     = evaluations.vec
        out$estim.precis    = estim.precis.vec
        }


    #   print( str(out) )

    if( FALSE )
        {
        time_elapsed = gettime() - time_start
        log.string( INFO, "Processed %d samples in %g seconds.  %g/sample.  [mean(time.elapsed)=%g].",
                                    n, time_elapsed, time_elapsed/n,  mean(out$time.elapsed,na.rm=TRUE) )
        }


    #print(mask)

    if( warn )
        {
        mask    = is.na(HVC[ ,1])  |  is.na(HVC[ ,3])
        if( any(mask)  )
            {
            log.string( WARN, "%d samples, out of %d, could not be mapped; HVC row set to NA.",
                                                sum(mask), n )
            }
        }

    return( out )
    }



#   this one use global variable p.InversionCoeffs
approxABfromxyY <- function( xyY, xyY.white, V.vector )
    {
    p   = 'spacesXYZ'
    if( ! requireNamespace( p, quietly=TRUE ) )
        {
        log.string( ERROR, "required package '%s' could not be loaded.", p )
        return(NULL)
        }

    V   = VfromY( xyY[3] )
    iV  = which.min( abs(V - V.vector) )

    coeffs  = p.InversionCoeffs[[iV]]
    vars    = colnames(coeffs)[1]

    if( grepl( "a,b", vars, fixed=TRUE ) )
        {
        #   polynomial in a and b
        Yp          = YfromV( V.vector[iV] )    # use Y in the closest V plane, which is where the model was computed
        XYZ         = spacesXYZ::XYZfromxyY( c(xyY[1],xyY[2],Yp) )
        XYZ.white   = spacesXYZ::XYZfromxyY( xyY.white )

        Lab = spacesXYZ::LabfromXYZ( XYZ, XYZ.white )
        a   = Lab[2]
        b   = Lab[3]
        a2  = a*a
        b2  = b*b
        term    = c( a, a2, a2*a, b, a*b, a2*b, b2, a*b2, b2*b )    # compare order here with polym() in makeInversionCoeffs()
        }
    else if( grepl( "xdelta,ydelta", vars, fixed=TRUE ) )
        {
        #   polynomial in xdelta and ydelta
        xd  = xyY[1] - xyY.white[1]
        yd  = xyY[2] - xyY.white[2]

        xd2 = xd*xd
        yd2 = yd*yd

        term    = c( xd, xd2, xd2*xd, yd, xd*yd, xd2*yd, yd2, xd*yd2, yd2*yd )    # compare order here with polym() in makeInversionCoeffs()
        }
    else
        log.string( FATAL, "vars='%s' is invalid.", vars )

    AB  = numeric(2)
    for( k in 1:2 )
        AB[k]   = sum( coeffs[k, ] * term )

    return( AB )
    }