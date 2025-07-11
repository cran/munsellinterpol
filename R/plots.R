


plotLociHC  <-  function( value=5, hue=seq(2.5,100,by=2.5), chroma='auto', coords='xy', main="Value %g/", est=FALSE, ... )
    {
    ok  = is.numeric(value)  &&  1<=length(value)  &&  all( 0<value  &  value<=10 )
    if( ! ok )
        {
        log_level( ERROR, "All Values must be numeric and in the interval (0,10]." )
        return(FALSE)
        }

    if( is.character(hue) )
        {
        hue = HueNumberFromString(hue)
        if( all( is.na(hue) ) )
            {
            log_level( "All character hues are invalid." )
            return(FALSE)
            }
        }

    ok  = is.numeric(hue)  &&  1<=length(hue)
    if( ! ok )
        {
        log_level( ERROR, "All Hues must be numeric." )
        return(FALSE)
        }

    hue = (hue %% 100)
    hue[ hue == 0 ] = 100

    ok  = (is.numeric(chroma)  &&  1<=length(chroma)  &&  all( 0 < chroma ))
    ok  = ok || (is.character(chroma) && chroma[1] == 'auto')
    if( ! ok )
        {
        log_level( ERROR, "All Chromas must be numeric and positive, or the string 'auto'." )
        return(FALSE)
        }

    if( ! (coords %in% c('xy','ab','AB') ) )
        {
        log_level( ERROR, "coords = '%s' is invalid.", coords )
        return(FALSE)
        }

    Value.vec   = sort( unique(munsellinterpol::Munsell2xy$V) )

    for( val in value )
        {
        main.val    = sprintf( main, val )

        if( is.character(chroma)  &&  chroma[1] == 'auto' )
            {
            idx     = which.min( abs(Value.vec - val) )
            vnear   = min( Value.vec[idx], 9 )
            mask    = (munsellinterpol::Munsell2xy$V == vnear)  &  (munsellinterpol::Munsell2xy$real)
            if( any(mask) )
                {
                chromamax   = max( munsellinterpol::Munsell2xy$C[mask] )
                chromavec  =  seq( 2, chromamax, by=2 )
                }
            else
                {
                log_level( ERROR, "Cannot determine maximum chroma for value=%g.", val )
                next
                }
            }
        else
            chromavec  = chroma

        if( coords == 'xy' )
            out = plotLociHC.xy( value, hue, chromavec, main.val, est, ...  )
        else if( coords == 'ab' )
            out = plotLociHC.ab( value, hue, chromavec, main.val, est, ...  )
        else if( coords == 'AB' )
            out = plotLociHC.AB( value, hue, chromavec, main.val, est, ...  )

        if( ! out ) break
        }

    return( invisible(out) )
    }



#   value       a *single* value this time
plotLociHC.xy  <-  function( value, hue, chroma, main, est, ... )
    {
    p   = 'spacesXYZ'
    if( ! requireNamespace( p, quietly=TRUE ) )
        {
        log_level( ERROR, "required package '%s' could not be loaded.", p )
        return(NULL)
        }

    df  = expand.grid( H=hue, C=chroma )

    HVC = as.matrix( cbind( df$H, value, df$C ) )

    #   add the illuminant at the end
    HVC = rbind( HVC, c(0,10,0) )

    xyY = MunsellToxyY( HVC, warn=FALSE, ... )$xyY

    xy_white    = xyY[ nrow(xyY), 1:2 ]

    xlim    = range( xyY[ ,1], na.rm=TRUE )
    ylim    = range( xyY[ ,2], na.rm=TRUE )



    plot( xlim, ylim, type='n', xlab='', ylab='', las=1, asp=1, lab=c(10, 10, 7), tcl=0, mgp=c(3,0.25,0) )
    title( xlab="x", line=1.25 )
    title( ylab="y", line=2.5 )
    grid( lty=1, col='gray70' )
    abline( h=0, v=0 )
    #   lines( c(0,1), c(1,0), lty=2 )

    #   draw the MacAdam limits
    tmp = sectionOptimals( YfromV(value) )
    if( ! is.null(tmp) )
        {
        xyY = spacesXYZ::xyYfromXYZ( tmp$XYZ )
        polygon( xyY[ ,1], xyY[ ,2], border='red', lwd=0.5 )
        }

    #   draw the inverted U
    xyz = p.xyz1931[  , 2:4 ]
    xyY = spacesXYZ::xyYfromXYZ( as.matrix(xyz) )
    polygon( xyY[ ,1], xyY[ ,2], border='blue', lwd=0.5 )


    #   put a small point at Illuminant C
    #xyY = MunsellToxyY( 'N 10/', warn=FALSE, ... )$xyY
    #xy_white = xyY[1:2]
    points( xy_white[1], xy_white[2], pch=19, cex=0.25 )

    #   draw radials
    chroma_min  = min( chroma, 1 )
    chroma_rad  = seq( chroma_min, max(chroma), by=0.25 )
    #   n   = length(chroma_rad)

    for( h in hue )
        {
        HVC = cbind( h, value, chroma_rad )

        xyY = MunsellToxyY( HVC, warn=FALSE, ...  )$xyY

        lines( xyY[ ,1], xyY[ ,2], lwd=0.5 )   #, border='black' )

        #   put a label at last point that is not NA
        #   print(  xyY[ ,1] )
        idx     = which( is.finite(xyY[ ,1]) )  #; print( is.finite(xyY[ ,1]) ) ; print(idx)
        m   = length(idx)
        xy      = xyY[idx[m],1:2]
        #offset  = xy - xy_white
        offset  = xy - xyY[idx[m-1],1:2]

        adj = 0.5 - 0.7 * sign(offset)
        text( xy[1], xy[2], HueStringFromNumber(h), adj=adj, cex=0.75, xpd=NA )   #, border='black' )
        }

    #   draw ovoids
    huevec  = seq( 0, 100, by=2.5/8 )
    #   n       = length(huevec)

    for( ch in chroma )
        {
        HVC = cbind( huevec, value, ch )

        xyY = MunsellToxyY( HVC, warn=FALSE, ... )$xyY

        lines( xyY[ ,1], xyY[ ,2], lwd=0.5 )   #, border='black' )
        }


    # label the hue lines
    #   n   = length(hue)
    #ch  = 1.05 * max(chroma)
    #HVC = cbind( hue, value, ch )
    #xyY = MunsellToxyY( HVC, warn=FALSE, ... )$xyY
    #text( xyY[ ,1], xyY[ ,2], HueStringFromNumber(hue), adj=c(0.5,0.5), cex=0.75 )   #, border='black' )


    #   mark the points in the lookup table
    mask    = (munsellinterpol::Munsell2xy$V == value)  &  (munsellinterpol::Munsell2xy$H %in% hue)  &  (munsellinterpol::Munsell2xy$C %in% chroma)

    if( any(mask) )
        {
        if( est )
            {
            data   = makeDataForPrediction( munsellinterpol::Munsell2xy, value, p.LookupList )  #; print( str(dfsub) )
            if( is.null(data) )    return(FALSE)

            #   print( dfsub )

            data   = addPredictions( data, warn=FALSE )

            #   for debugging
            #   print( data[ nrow(data), ] )

            segments( data$x, data$y, data$x.pred, data$y.pred, col='red' )
            points( data$x.pred, data$y.pred, pch=19, cex=0.4, col='red' )
            }

        dfsub   = munsellinterpol::Munsell2xy[ mask,  ]
        pch     = ifelse( dfsub$real, 19, 1 )
        points( dfsub$x, dfsub$y, pch=pch, cex=0.5, lwd=0.25 )   #, border='black' )
        }


    mess = sprintf( "%s \n[Chromas = %s]", main, paste0(chroma,collapse=', ') )
    #mess    = main
    title( main=mess )

    return(TRUE)
    }



#   value   a *single* value
#   hue     can be a vector
#   chroma  can be a vector

plotLociHC.ab  <-  function( value, hue, chroma, main, est, ...  )
    {
    p   = 'spacesXYZ'
    if( ! requireNamespace( p, quietly=TRUE ) )
        {
        log_level( ERROR, "required package '%s' could not be loaded.", p )
        return(NULL)
        }

    # par( omi=rep(0.5,4) )

    XYZ.C   = spacesXYZ::XYZfromxyY( c( p.xyC['NBS', ], 100 ) )

    df  = expand.grid( H=hue, C=chroma )

    HVC = as.matrix( cbind( df$H, value, df$C ) )

    #   add the illuminant at the end
    HVC = rbind( HVC, c(0,10,0) )

    xyY = MunsellToxyY( HVC, warn=FALSE, ... )$xyY

    #   xy_white    = xyY[ nrow(xyY), 1:2 ]

    #   convert to XYZ
    XYZ     = spacesXYZ::XYZfromxyY( xyY )
    XYZ.C   = XYZ[ nrow(XYZ), ]

    Lab     = spacesXYZ::LabfromXYZ( XYZ, XYZ.C )

    xlim    = range( Lab[ ,2], na.rm=TRUE )
    ylim    = range( Lab[ ,3], na.rm=TRUE )


    plot( xlim, ylim, type='n', xlab='', ylab='', las=1, asp=1, lab=c(10, 10, 7),  tcl=0, mgp=c(3,0.25,0)  )
    title( xlab="a", line=1.25 )
    title( ylab="b", line=2.5 )
    grid( lty=1, col='gray70' )
    abline( h=0, v=0 )

    #   draw the MacAdam limits
    tmp = sectionOptimals( YfromV(value) )
    if( ! is.null(tmp) )
        {
        Lab = spacesXYZ::LabfromXYZ( tmp$XYZ, XYZ.C )
        polygon( Lab[ ,2], Lab[ ,3], border='red', lwd=0.5 )
        }

    #   draw radials
    # cat( "drawing radials...\n", file=stderr() )

    chromavec  = seq( 1, max(chroma), by=0.25 )
    #   n   = length(chroma)

    for( h in hue )
        {
        HVC = cbind( h, value, chromavec  )

        xyY = MunsellToxyY( HVC, warn=FALSE, ...  )$xyY
        XYZ = spacesXYZ::XYZfromxyY( xyY )
        Lab = spacesXYZ::LabfromXYZ( XYZ, XYZ.C )

        lines( Lab[ ,2], Lab[ ,3], lwd=0.5 )   #, border='black' )

        #   put a label at the end
        idx     = which( is.finite(Lab[ ,2]) )  #; print( is.finite(xyY[ ,1]) ) ; print(idx)
        m       = length(idx)
        ab      = Lab[idx[m],2:3]
        offset  = ab - Lab[idx[m-1],2:3]

        adj = 0.5 - 0.7 * sign(offset)
        text( ab[1], ab[2], HueStringFromNumber(h), adj=adj, cex=0.75, xpd=NA )   #, border='black' )

        #  if( any(is.na(XYZ)) )   print( xyY )
        }



    #   draw ovoids
    # cat( "drawing ovoids...\n", file=stderr() )
    huevec  = seq( 0, 100, by=2.5/8 )
    #   n       = length(huevec)
    for( ch in chroma )
        {
        HVC = cbind( huevec, value, ch )

        xyY = MunsellToxyY( HVC, warn=FALSE, ... )$xyY
        XYZ = spacesXYZ::XYZfromxyY( xyY )
        Lab = spacesXYZ::LabfromXYZ( XYZ, XYZ.C )

        lines( Lab[ ,2], Lab[ ,3], lwd=0.5  )   #, border='black' )
        }


    #   mark the points in the lookup table

    data   = makeDataForPrediction( munsellinterpol::Munsell2xy, value, p.LookupList )  #; print( str(dfsub) )

    if( ! is.null(data) )
        {
        mask    =  (data$H %in% hue)  &  (data$C %in% chroma)

        if( any(mask)  &&  est )
            {
            data   = addPredictions( data, warn=FALSE )

            dfsub   = data[ mask, , drop=FALSE ]

            segments( dfsub$a, dfsub$b, dfsub$a.pred, dfsub$b.pred, col='red' )
            points( dfsub$a.pred, dfsub$b.pred, pch=19, cex=0.4, col='red' )
            }

        pch     = ifelse( data$real, 19, 1 )
        points( data$a, data$b, pch=pch, cex=0.5, lwd=0.25 )   #, border='black' )
        }


    mess = sprintf( "%s \n[Chromas=%s]", main, paste0( chroma, collapse=', ' ) )
    title( main=mess )

    return( invisible(TRUE) )
    }



plotLociHC.AB  <-  function( value, hue, chroma, main, est, ...    )
    {
    chroma.max  = max(chroma)

    xlim    = c( -chroma.max, chroma.max )
    ylim    = c( -chroma.max, chroma.max )

    plot( xlim, ylim, type='n', xlab='A', ylab='B', las=1, asp=1, lab=c(10, 10, 7) )


    #   draw radials
    for( theta in (hue * pi/50) )
        lines( chroma.max * c(0,cos(theta)), chroma.max * c(0,sin(theta)), col="#eeeeee" )   #, border='black' )

    #   draw circles
    theta   = seq( 0, 2*pi, len=101 )
    circle  = matrix( c(cos(theta),sin(theta)), length(theta), 2 )
    for( ch in chroma )
        lines( ch * circle[ ,1], ch * circle[ ,2], col="#eeeeee" )

    grid( lty=1 )
    abline( h=0, v=0 )


    data   = makeDataForPrediction( munsellinterpol::Munsell2xy, value, p.LookupList )  #; print( str(dfsub) )

    if( ! is.null(data) )
        {
        mask    =  (data$H %in% hue)  &  (data$C %in% chroma)

        if( est  &&  any(mask)  )
            {
            #    plot inversion estimates
            data   = addPredictions( data, warn=FALSE )

            dfsub   = data[ mask, , drop=FALSE ]

            segments( dfsub$A, dfsub$B, dfsub$A.pred, dfsub$B.pred, col='red' )
            points( dfsub$A.pred, dfsub$B.pred, pch=19, cex=0.4, col='red' )
            }

        pch     = ifelse( data$real, 19, 1 )
        points( data$A, data$B, pch=pch, cex=0.5, lwd=0.25 )
        }


    if( 0 )
    {
    #   plot the "crude" prediction
    dfreal   = dfsub[ dfsub$real, ]

    A.crude = dfsub$a / 5.5
    B.crude = dfsub$b / 5.5
    segments( dfreal$A, dfreal$B, A.crude, B.crude, col='blue' )
    points( A.crude, B.crude, pch=19, cex=0.4, col='blue' )

    #   plot the real points
    points( dfreal$A, dfreal$B, pch=19, cex=0.4 )

    #   plot the non-real points
    dfnonreal   = dfsub[ ! dfsub$real, ]
    points( dfnonreal$A, dfnonreal$B, pch=1, cex=0.4, col='cyan' )
    }


    mess = sprintf( "%s \n[Chromas=%s]", main, paste( chroma, collapse=' ' ) )
    title( main=mess )

    return(TRUE)
    }


plotPatchesH  <-  function( hue, space='sRGB', adapt='Bradford', background='gray50',
                                main="Hue %s  (H=%g)      [%s   adapt=%s]", value=NULL, chroma=NULL, ...  )
    {
    if( is.character(hue) )
        {
        hue = HueNumberFromString(hue)
        if( any( is.na(hue) ) )
            {
            log_level( "One or more of the character hues is invalid." )
            return( invisible(FALSE) )
            }
        }

    ok  = is.numeric(hue)  &&  1<=length(hue)
    if( ! ok )
        {
        log_level( ERROR, "All Hues must be numeric." )
        return( invisible(FALSE) )
        }

    hue = (hue %% 100)
    hue[ hue == 0 ] = 100

    if( ! is.null(value) )
        {
        #   check value
        ok  = is.numeric(value)  &&  0<length(value)  &&  0<=value[1]  &&  all( 0 < diff(value) )  &&  value[length(value)]<=10
        if( ! ok )
            {
            log_level( ERROR, "argument 'value' is invalid.  It must be strictly increasing with all values in [0,10]." )
            return( invisible(FALSE) )
            }
        }

    if( ! is.null(chroma) )
        {
        #   check chroma
        ok  = is.numeric(chroma)  &&  0<length(chroma)  &&  0<=chroma[1]  &&  all( 0 < diff(chroma) )
        if( ! ok )
            {
            log_level( ERROR, "argument 'chroma' is invalid.  It must be non-negative and strictly increasing." )
            return( invisible(FALSE) )
            }
        }

    for( h in hue )
        {
        main.single = sprintf( main, HueStringFromNumber(h), h, space, adapt )

        out = plotPatchesH.single( h, space=space, adapt=adapt, background=background, main=main.single,
                                        valuevec=value, chromavec=chroma, ... )

        if( ! out ) break
        }

    return( invisible(out) )
    }




#   hue is a single number
plotPatchesH.single <- function( hue, space='sRGB', adapt='Bradford', background='gray50', main='',
                                        valuevec=NULL, chromavec=NULL, ... )
    {
    #   we might have to compute HVC twice, once from real.dat and once as a grid
    makegrid    = ! is.null(valuevec)  ||  ! is.null(chromavec)

    if( is.null(valuevec) )
        valuevec = 0:10

    if( is.null(chromavec) )
        {
        #   find the closest hue page
        huevec  = sort( unique( munsellinterpol::Munsell2xy$H ) )   #; print( huevec )

        delta   = abs(hue - huevec)
        delta   = pmin( delta, 100 - delta )
        hue.near    = huevec[ which.min( delta ) ]  #; print( hue.near )

        #   find all the patches on the hue.near page, using samples in the file real.dat
        #   only look at real patches, and Value>=1
        mask    = munsellinterpol::Munsell2xy$H == hue.near  &  munsellinterpol::Munsell2xy$real  &  1 <= munsellinterpol::Munsell2xy$V
        dfsub   = munsellinterpol::Munsell2xy[ mask, ]

        HVC     = cbind( hue, dfsub$V, dfsub$C )

        if( ! makegrid )
            #   the HVC we are making now is the final one, but it currently has no neutrals, so add them
            HVC = rbind( cbind(0,valuevec,0), HVC )

        #   find the largest Chroma in the hue.near page, using samples in  the file real.dat
        chromamax   = max( HVC[ ,3] )

        chromavec   = seq( 0, chromamax, by=2 )
        }

    if( makegrid )
        {
        #   ignore previous HVC and create a new one
        #   the previous one was only needed to get chromavec
        df  = expand.grid( valuevec, chromavec )

        HVC = cbind( hue, df[ ,1], df[ ,2] )
        }

    colnames(HVC)  = c('H','V','C')

    #   the next call might create WARNs, so suppress them
    tmp     = MunsellToRGB( HVC, space=space, adapt=adapt, warn=FALSE, ... )

    #cat( "Bad xyY's = ", sum( is.na(tmp$xyY[ ,1]) ), '\n' )
    #print( tmp )

    xyY     = tmp$xyY          # this is adapted to illuminant C

    #   only keep ones that were successfully transformed, and then within MacAdamLimits
    tol     = 0.004     # adjusted for "N 10/", which is problematic, it's right on the boundary
    mask    = is.finite( xyY[ ,1] )   &   (IsWithinMacAdamLimits( xyY, 'C' )$delta <= tol)     #;   print( sum(mask) )

    # the sample N 10/ is problematic, it's right on the boundary
    #idx = which( HVC[ ,2]==10  &  HVC[ ,3]==0 )
    #if( 0 < length(idx) )   mask[ idx[1] ]  = TRUE

    #   throw away the ones outside the limits
    dfsub       = tmp[ mask, , drop=FALSE ]
    dfsub$HVC   = HVC[ mask, , drop=FALSE ]

    # print( dfsub )
    # cat( "samples inside MacAdamLimits = ", nrow(dfsub), '\n' )

    dfsub$color = grDevices::rgb( round(dfsub$RGB), maxColorValue=255 )    #    print( str(dfsub) )

    #   setup the plot
    xlim    = c( 0.5, length(chromavec)+0.5 )
    ylim    = c( 0.5, length(valuevec)+0.5 )

    bg.prev = par( bg=background )
    par( mgp=c(0, 0.5, 0) )
    plot.new()
    plot.window( xlim, ylim )


    #   set up margins
    margin  = 0.1

    size    = grDevices::dev.size('in')

    scale.x = size[1] / (xlim[2] - xlim[1])
    scale.y = size[2] / (ylim[2] - ylim[1])


    if( scale.y < scale.x )
        {
        margin.x    = margin * (scale.y / scale.x)
        margin.y    = margin
        }
    else
        {
        margin.x    = margin
        margin.y    = margin * (scale.x / scale.y)
        }

    #   draw axes
    lines( c(xlim[1],xlim[2]) - 0, c(ylim[1],ylim[1]) - 0, xpd=NA )   # X axis
    lines( c(xlim[1],xlim[1]) - 0, c(ylim[1],ylim[2]) - 0, xpd=NA )   # Y axis

    title( xlab='Chroma', line=0.5 )
    title( ylab='Value', line=0.5 )

    title( main=main, line=0 )

    #   x labels
    text( 1:length(chromavec), rep(ylim[1],length(chromavec)), sprintf( "/%g", chromavec ), adj=c(0.5,1.5), cex=0.8, xpd=NA )

    #   y labels  =  Value labels
    ht  = graphics::strheight('0')
    text( rep(xlim[1] - ht/5, length(valuevec)), 1:length(valuevec), sprintf( "%g/", valuevec ), adj=c(1,0.5), cex=0.8, xpd=NA )

    if( nrow(dfsub) == 0 )
        {
        #   nothing to draw
        par( bg=bg.prev )   # restore previous background
        return(FALSE)
        }

    #   make small dictionaries for Value and Chroma
    xfromchroma = 1:length(chromavec)
    names(xfromchroma)  = as.character(chromavec)

    yfromvalue  = 1:length(valuevec)
    names(yfromvalue)  = as.character(valuevec)

    if( FALSE &&  chromavec[1] == 0 )
        {
        #   draw neutrals
        HVC = cbind( 0, valuevec, 0 )
        dfneutral   = MunsellToRGB( HVC, space=space, adapt=adapt, ... )
        dfneutral$HVC   = HVC
        dfneutral$color = grDevices::rgb( round(dfneutral$RGB), maxColorValue=255 )   #    print( dfneutral )

        n       = nrow( dfneutral )
        left    = rep( 0.5 + margin.x, n )
        right   = rep( 1.5 - margin.x, n )
        y       = yfromvalue[ as.character(dfneutral$HVC[ ,2]) ]
        bottom  = y - 1/2 + margin.y
        top     = y + 1/2 - margin.y
        rect( left, bottom, right, top, col=dfneutral$color, border=NA )
        }


    #   draw the HVCs
    x   = xfromchroma[ as.character(dfsub$HVC[ ,3]) ]
    y   = yfromvalue[ as.character(dfsub$HVC[ ,2]) ]

    dfsub$left    = x - 1/2 + margin.x
    dfsub$right   = x + 1/2 - margin.x
    dfsub$bottom  = y - 1/2 + margin.y
    dfsub$top     = y + 1/2 - margin.y
    color   = ifelse( dfsub$OutOfGamut, NA, dfsub$color )
    border  = ifelse( dfsub$OutOfGamut, 'gray30', NA )
    rect( dfsub$left, dfsub$bottom, dfsub$right, dfsub$top, col=color, border=border )

    #   draw RGB labels inside the out-of-RGB-gamuts
    dfout   = dfsub[ dfsub$OutOfGamut, , drop=FALSE ]
    if( 0 < nrow(dfout) )
        {
        RGB     = round( dfout$RGB )
        label   = sprintf( "R %3d\nG %3d\nB %3d", RGB[ ,1], RGB[ ,2], RGB[ ,3] )

        w   = graphics::strwidth( label[1], family="mono" )
        h   = graphics::strheight( label[1], family="mono" )
        w0  = graphics::strwidth( "0", family="mono" )

        cex = min( (1 - 2*margin.x - w0) / w, (1 - 2*margin.y - w0) / h )

        text( (dfout$left+dfout$right)/2, (dfout$top+dfout$bottom)/2, label, cex=cex, adj=c(0.5,0.5), family="mono" )
        }

    par( bg=bg.prev )   # restore previous background

    return(TRUE)
    }




plotLookupDF <- function( Munsell2xy )
    {
    Value.vec   = sort( unique(Munsell2xy$V) )

    huelim    = c(0,100)
    chromalim    = c(0,50)

    hueseq      = seq( 0, 100, by=2.5 )
    chromaseq   = seq( 0, 50, by=2 )
    for( v in Value.vec )
        {
        plot.new()
        plot.window( huelim, chromalim ) #,  mgp=c(2,2,2)  )
        #   par( mai=c( 0.5, 0.5, 0.5, 0 ) )
        #par( mgp=c(1/2,3,1/2) )
        title( xlab="Hue", ylab="Chroma", main=sprintf( "Value %g", v ) )

        segments(  0, chromaseq, 100, chromaseq, col='gray' )
        segments(  hueseq, 0, hueseq, 50, col='gray' )
        rect( huelim[1], chromalim[1], huelim[2], chromalim[2], border='black' )

        text( -1, chromaseq, as.character(chromaseq), adj=c(1,1/2), cex=0.75 )
        text( hueseq, -1, HueStringFromNumber(hueseq), adj=c(1,1/2), cex=0.75, srt=90, xpd=NA )

        dfsub   =   Munsell2xy[ Munsell2xy$V==v,  ]

        pch = ifelse( dfsub$real, 19, 1 )
        points( dfsub$H, dfsub$C, pch=pch )
        #   break
        }

    return(TRUE)
    }

plotLookupList <- function( LookupList, Munsell2xy, value=NULL )
    {
    Value.vec   = attr( LookupList, 'V.vector' )
    Hue.vec     = attr( LookupList, 'H.vector' )
    Chroma.vec  = attr( LookupList, 'C.vector' )

    huelim    = c(0,100)
    chromalim    = c(0,50)

    hueseq      = seq( 0, 100, by=2.5 )
    chromaseq   = seq( 0, 50, by=2 )

    if( is.null(value) )
        vvec    =  1:length(LookupList)
    else
        {
        vvec    = match( value, Value.vec )
        vvec    = vvec[ is.finite(vvec) ]
        }

    for( iV in vvec )
        {
        par( mai=c( 0.5, 0.5, 0.75, 0.25 ) )
        #par( mgp=c(1/2,3,1/2) )
        plot.new()
        plot.window( huelim, chromalim ) #,  mgp=c(2,2,2)  )

        # plot.default( huelim, chromalim, las=1, xlab='', ylab='', type='n', lab=c(10,8,7),  tcl=0, mgp=c(3, 0.25, 0) )

        title( xlab="Hue", line=0.25  )
        title( ylab="Chroma", line=0  )

        value   = Value.vec[iV]
        title(  main=sprintf( "Value %g", value ) )

        segments(  0, chromaseq, 100, chromaseq, col='gray' )
        segments(  hueseq, 0, hueseq, 50, col='gray' )
        rect( huelim[1], chromalim[1], huelim[2], chromalim[2], border='black' )

        text( -1, chromaseq, as.character(chromaseq), adj=c(1,1/2), cex=0.75 )
        text( 100+1, chromaseq, as.character(chromaseq), adj=c(0,1/2), cex=0.75 )

        hueskip = seq(0,100,by=5)
        text( hueskip, -1, hueskip, adj=c(0.5,1), cex=0.75 )
        text( hueseq, chromalim[2]+0.5, HueStringFromNumber(hueseq), adj=c(0,1/2), cex=0.75, srt=90, xpd=NA )

        dfsub   = Munsell2xy[ Munsell2xy$V==value,  ]
        pch = ifelse( dfsub$real, 19, 1 )
        points( dfsub$H, dfsub$C, pch=pch )

        dfsub   = dfsub[ dfsub$H==100, ]
        pch = ifelse( dfsub$real, 19, 1 )
        points( rep(0,length(dfsub$C)), dfsub$C, pch=pch )

        idx = which( ! is.na(LookupList[[iV]]$x), arr.ind=TRUE )  #; print( idx )
        hue     = Hue.vec[ idx[ ,1] ]
        chroma  = Chroma.vec[ idx[ ,2] ]
        mask    = 0 <= hue  &  hue <= 100
        points( hue[mask], chroma[mask], pch=4 )
        #   break

        #dfsub   = subset( Munsell2xy, V==v )

        #pch = ifelse( dfsub$real, 19d, 1 )
        #points( dfsub$H, dfsub$C, pch=pch )
        #   break
        }

    return(TRUE)
    }

#   x   mask is TRUE
#       solid disk, real
#   o   open circle, not real
plotMask3D <- function( Mask3D, Munsell2xy )
    {
    Hue.vec     = as.numeric( dimnames(Mask3D)[[1]] )
    Value.vec   = as.numeric( dimnames(Mask3D)[[2]] )
    Chroma.vec  = as.numeric( dimnames(Mask3D)[[3]] )

    huelim    = c(0,100)
    chromalim    = c(0,50)

    hueseq      = seq( 0, 100, by=2.5 )
    chromaseq   = seq( 0, 50, by=2 )


    for( iV in 1:length(Value.vec) )
        {
        plot.new()
        plot.window( huelim, chromalim ) #,  mgp=c(2,2,2)  )
        #   par( mai=c( 0.5, 0.5, 0.5, 0 ) )
        #par( mgp=c(1/2,3,1/2) )

        value   = Value.vec[iV]
        title( xlab="Hue", ylab="Chroma", main=sprintf( "Value %g", value ) )

        segments(  0, chromaseq, 100, chromaseq, col='gray' )
        segments(  hueseq, 0, hueseq, 50, col='gray' )
        rect( huelim[1], chromalim[1], huelim[2], chromalim[2], border='black' )

        text( -1, chromaseq, as.character(chromaseq), adj=c(1,1/2), cex=0.75 )
        text( hueseq, -1, HueStringFromNumber(hueseq), adj=c(1,1/2), cex=0.75, srt=90, xpd=NA )


        dfsub   = Munsell2xy[ Munsell2xy$V==value,  ]

        pch = ifelse( dfsub$real, 19, 1 )
        points( dfsub$H, dfsub$C, pch=pch )

        idx = which( Mask3D[ , iV, ], arr.ind=TRUE )
        H   = Hue.vec[ idx[ ,1] ]
        C   = Chroma.vec[ idx[ ,2] ]
        points( H, C, pch=4 )
        }

    return(TRUE)
    }


