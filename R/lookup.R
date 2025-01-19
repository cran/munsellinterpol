
#   p.LookupList    = NULL

makeLookupList <- function( Munsell2xy, xyC=c(0.3101,0.3163), kfactor=c(0.7,0.5) )
    {   
    time_start  = gettime()
    
    Hue.vector      = seq(2.5,100,by=2.5)   #  H is always a multiple of 2.5. later this vector will change
            
    Value.vector    = sort( unique(Munsell2xy$V) )
    Value.vector    = c( 0, Value.vector )  # put 0 at beginning, so we can interpolate between 0 and 0.2
    #   Value.vector    = c(0,0.2,0.4,0.6,0.8,1,2,3,4,5,6,7,8,9,10)
    
    #   put 0 at beginning of chroma.vector,
    #   because we will put the neutral x,y there, for Illuminant C
    Chroma.vector   = seq(0,50,by=2) 
    
    #   allocate 4D array, full of NAs
    Lookup4D = array( NA_real_, dim=c(length(Hue.vector),length(Value.vector),length(Chroma.vector),2) )
    
    dimnames(Lookup4D)[[1]] = as.character( Hue.vector )    
    dimnames(Lookup4D)[[2]] = as.character( Value.vector )
    dimnames(Lookup4D)[[3]] = as.character( Chroma.vector )
    dimnames(Lookup4D)[[4]] = c('x','y')
    
    #   filter out very dark samples that are not real
    #mask    = Munsell2xy$real  |  (1 <= Munsell2xy$V )
    #dfsub   = Munsell2xy[ mask, ]
    
    #   copy Munsell2xy to Lookup4D[].  Some of these may be overwritten later   
    for( i in 1:nrow(Munsell2xy) )
        {
        row = Munsell2xy[i, ]
        
        iH  = match( row$H, Hue.vector )
        iV  = match( row$V, Value.vector )
        iC  = match( row$C, Chroma.vector )        
        
        Lookup4D[ iH, iV, iC,   ] = c( row$x, row$y )
        }    


    
    Lookup4D = extrapGAM( Lookup4D, Munsell2xy, kfactor=kfactor  )

    
    #print( str(Lookup4D) ) ;  print( object.size(Lookup4D) )
    #   return(FALSE)
    
    #   count the total number of non-NA points in Lookup4D
    gridcount   = integer( length(Value.vector) )
    names(gridcount)    = as.character(Value.vector)
    
    for( iV in 2:length(Value.vector) )
        gridcount[iV]   = sum( is.finite(Lookup4D[ ,iV, ,'x',drop=FALSE]) )
    
    gridcount   = c( gridcount, sum(gridcount) )
    names(gridcount)[ length(gridcount) ]   = 'sum'
    cat( "Total gridcount:\n" )
    print( gridcount )

    
    
    #   now copy from array Lookup4D[] to a list with an item for each Munsell Value
        
    out = vector( length(Value.vector), mode='list' )
    
    names(out) = as.character(Value.vector)

     #   there are are only 40 hues, but we make 43 for easier wraparound (periodicity)
    Hue.vector      = seq(-2.5,102.5,by=2.5) 
        
    attr( out, "V.vector" ) = Value.vector
    attr( out, "H.vector" ) = Hue.vector
    attr( out, "C.vector" ) = Chroma.vector
    
    #   initialize with arrays of NAs
   
    dummy   = array( NA_real_, dim=c( length(Hue.vector), length(Chroma.vector)  ) )
    rownames(dummy) = HueStringFromNumber( Hue.vector )
    colnames(dummy) = as.character( Chroma.vector )
    
    for( iV in 1:length(out) )
        {
        out[[iV]] = list()  # list( x=dummy, y=dummy )    
        
        dummy[ 3:42, ]  = Lookup4D[ ,iV, , 'x' ]
        out[[iV]]$x     = dummy        
        
        dummy[ 3:42, ]  = Lookup4D[ ,iV, , 'y' ]
        out[[iV]]$y     = dummy        
        }
    
        
    #   copy data from V=0.2 to V=0
    out[[1]]$x = out[[2]]$x
    out[[1]]$y = out[[2]]$y
    
    #   trim away high Chroma, to reduce memory
    for( iV in 1:length(out) )
        {
        x   = out[[iV]]$x  #; print( str(x) )
        
        cmax    = NA        
        for( iC in 2:ncol(x) )
            {
            if( all( is.na(x[ ,iC]) ) )
                {
                cmax = iC - 1
                break
                }
            }
            
        if( is.na(cmax) )
            {
            log_level( INFO, "At V=%g, there is data up to the maximum C=%g.",
                            Value.vector[iV], Chroma.vector[length(Chroma.vector)] )
            attr( out[[iV]], "C.vector" ) = Chroma.vector            
            }
        else
            {
            log_level( INFO, "At V=%g, there is data up to C=%g.",
                                Value.vector[iV], Chroma.vector[cmax] )
                
            out[[iV]]$x = out[[iV]]$x[ , 1:cmax ]
            out[[iV]]$y = out[[iV]]$y[ , 1:cmax ]
            
            attr( out[[iV]], "C.vector" ) = Chroma.vector[ 1:cmax ]            
            }
            
            
        #   duplicate data for Hue wraparound
        out[[iV]]$x[1, ] = out[[iV]]$x[41, ]
        out[[iV]]$y[1, ] = out[[iV]]$y[41, ]
            
        out[[iV]]$x[2, ] = out[[iV]]$x[42, ]
        out[[iV]]$y[2, ] = out[[iV]]$y[42, ]    

        out[[iV]]$x[43, ] = out[[iV]]$x[3, ]
        out[[iV]]$y[43, ] = out[[iV]]$y[3, ]    
        
        #   replace NAs by 0s
        #out[[iV]]$x[ is.na(out[[iV]]$x) ] = 0
        #out[[iV]]$y[ is.na(out[[iV]]$y) ] = 0      
        }
        
        
    # fill Chroma=0 with xyC
    for( k in 1:length(out) )
        {
        out[[k]]$x[ ,1] = xyC[1]
        out[[k]]$y[ ,1] = xyC[2]
        }        
    
    time_elapsed  = gettime()  -  time_start
    
    log_level( INFO, "Made list in %g seconds.", time_elapsed )
    
    return( invisible(out) )
    }

    
extrapGAM  <-  function( Lookup4D, Munsell2xy, kfactor=c(0.9,0.5)  )
    {
    if( ! requireNamespace( 'mgcv', quietly=TRUE ) )        return( NULL )
    
    if( ! requireNamespace( 'spacesXYZ', quietly=TRUE ) )   return(NULL)
        
    
    time_start  = gettime()
    
    Hue.vector      = as.numeric( dimnames(Lookup4D)[[1]] )
    Value.vector    = as.numeric( dimnames(Lookup4D)[[2]] )    
    Chroma.vector   = as.numeric( dimnames(Lookup4D)[[3]] )
    
    
    out = Lookup4D
    
    Mask3D  = makeRequiredMask3D( Munsell2xy )    
    
    samples   = 0
    
    log_level( INFO, "Processing colors with V < 1..."  )
    
    #   for value < 1 make a special data.frame, and make only 1 model from it
    #   use a subset of points for V<1, plus all points for V==1
    mask    = Munsell2xy$V < 1  
    mask    = mask &  Munsell2xy$real  #  |  (75 <= Munsell2xy$H  &  Munsell2xy$H <= 100) )
    mask    = mask | (Munsell2xy$V == 1)

    data.dark   = Munsell2xy[ mask, ]    

    #   add rectangular coordinates A,B
    tmp = ABfromHC( data.dark$H, data.dark$C )
    data.dark   = cbind( data.dark, A = tmp$A, B = tmp$B )
    
    #   add ab coordinates
    
    white.C = spacesXYZ::XYZfromxyY( c( p.xyC['NBS',],100) )
    Y       = YfromV( data.dark$V )
    XYZ     = spacesXYZ::XYZfromxyY( cbind( data.dark$x, data.dark$y, Y ) )
    Lab     = spacesXYZ::LabfromXYZ( XYZ, white.C )
    
    data.dark$a = Lab[ ,2]
    data.dark$b = Lab[ ,3]
    
    #Luv     = xyz2luv( XYZ, white.C )
    #data.dark$u = Luv[ ,2]
    #data.dark$v = Luv[ ,3]
    
    
    #   print( str(data.dark) )
    if( TRUE )
    {
    k   = round( kfactor[1] * nrow(data.dark) )
    
    form    = sprintf( "x ~ s( A, B, V,  bs='tp', k=%d )", k )
    log_level( INFO, "Computing model '%s'...", form ) ; flush.console();
    gam.x = mgcv::gam( formula(form), data=data.dark )  

    form    = sprintf( "y ~ s( A, B, V,  bs='tp', k=%d )", k )
    log_level( INFO, "Computing model '%s'...", form ) ; flush.console();    
    gam.y = mgcv::gam( formula(form), data=data.dark )    
    
    log_level( INFO, "Very Dark"  )
    log_level( INFO, "x residual range = [%g,%g]", range(residuals(gam.x))[1],  range(residuals(gam.x))[2] )
    log_level( INFO, "y residual range = [%g,%g]", range(residuals(gam.y))[1],  range(residuals(gam.y))[2] )
    }
    
    
    
    for( value in c(0.2,0.4,0.6,0.8) )
        {
        iV  = match( value, Value.vector )
        
        if( FALSE )
        {
        dfv = Munsell2xy[ Munsell2xy$V==value, ]

        #   make mask of the ones to keep for modeling
        mask.model  =   dfv$real    #   |  (75 <= dfv$H  &  dfv$H <= 100)
        
        dfmodel     = dfv[ mask.model, ]
        
        #   add rectangular coordinates A,B
        tmp = ABfromHC( dfmodel$H, dfmodel$C )
        dfmodel  = cbind( dfmodel, A=tmp$A, B=tmp$B )
        
        k   = round( kfactor[1] * nrow(dfmodel) )
        
        form    = sprintf( "x ~ s( H, C,  bs='tp', k=%d )", k )
        log_level( INFO, "value=%g.  Computing model '%s'...", value, form ) ; flush.console();
        gam.x = mgcv::gam( formula(form), data=dfmodel )  

        form    = sprintf( "y ~ s( H, C,  bs='tp', k=%d )", k )
        log_level( INFO, "value=%g.  Computing model '%s'...", value, form ) ; flush.console();    
        gam.y = mgcv::gam( formula(form), data=dfmodel )    

        log_level( INFO, "x residual range = [%g,%g]", range(residuals(gam.x))[1],  range(residuals(gam.x))[2] )
        log_level( INFO, "y residual range = [%g,%g]", range(residuals(gam.y))[1],  range(residuals(gam.y))[2] )
        }
        
        
        #log_level( INFO, "For value=%g, found %d real samples.", value, nrow(dfreal) )        
        
        #   make 2D mask of the ones used for modeling
        mask2D  = matrix( FALSE, length(Hue.vector), length(Chroma.vector) )
        
        for( i in 1:nrow(data.dark) )
            {
            if( data.dark$V[i] == value )
                {
                iH  = match( data.dark$H[i], Hue.vector )
                iC  = match( data.dark$C[i], Chroma.vector )
                mask2D[iH,iC] = TRUE
                }
            }
        
        #   make a mask of iH,iC points to be extrapolated.  Preserve the ones used for modeling.
        #   preserve the entries used
        mask_extrap = Mask3D[ , iV, ]  &  (! mask2D)    # mask_extrap is 2D
        
        for( iC in 2:length(Chroma.vector) )    # start at 2 to ensure that neutral is NOT extrapolated
            {
            mask    = mask_extrap[ , iC]            # mask is 1D horizontal slice of constant C
            
            if( ! any(mask) )
                {
                #   log_level( INFO, "For value=%g, and iC=%d, no extrapolation needed.", value, iC )        
                next    # no extrapolation needed here
                }
                
            newdata = data.frame( H=Hue.vector[mask], V=value, C=Chroma.vector[iC] )
            
            tmp = ABfromHC( newdata$H, newdata$C )            
            newdata = cbind( newdata, A = tmp$A, B = tmp$B )
            
            #L   = Lightness_from_linear( YfromV( newdata$V ) / 100 )
            x   = predict( gam.x, newdata=newdata )
            y   = predict( gam.y, newdata=newdata )
            
            #Luv = cbind( L, u, v )
            #xyY = XYZ2xyY( luv2xyz( Luv, white.C ) )
            xyY = cbind( x, y, YfromV( newdata$V ) / 100 )
            
            mask.na = is.na( x )
            if( any(mask.na) )
                {
                log_level( WARN, "For value=%g and chroma=%g,  %d of the extrapolated xy's are NA.",
                                        value,  Chroma.vector[iC], sum(mask.na) )
                df  = newdata[mask.na, , drop=FALSE]
                #   df$Luv  = Luv[mask.na, , drop=FALSE]
                df$xyY  = xyY[mask.na, , drop=FALSE]
                print(df)
                }
            
            
            out[mask,iV,iC,'x']  =  x   # predict( gam.x, newdata=newdata )
            out[mask,iV,iC,'y']  =  y   # predict( gam.b, newdata=newdata )
            
            samples = samples + sum(mask)
            }   
        }
    
    log_level( INFO, "Extrapolated %d very dark samples in %g seconds.", samples, gettime() - time_start )  

    #   for V >= 1 do one value at a time
    log_level( INFO, "Processing colors with V >= 1..." )
    
    for( value in 1:10 )
        {
        iV  = match( value, Value.vector )

        #   for value >=1 Munsell2xy has good data
        dfsub   = Munsell2xy[ Munsell2xy$V==value,  ]  #; print(value) ; print( str(dfsub) )
        
        #   add rectangular coordinates A,B
        tmp = ABfromHC( dfsub$H, dfsub$C )                 
        dfsub   = cbind( dfsub, A = tmp$A, B = tmp$B ) 
            
        k   = round( kfactor[2] * nrow(dfsub) )  # - 10
        
        mess    = sprintf( "x ~ s( A, B, bs='tp', k=%d )", k )        
        gam.x   = mgcv::gam( formula(mess), data=dfsub )
        mess    = sprintf( "y ~ s( A, B, bs='tp', k=%d )", k )          
        gam.y   = mgcv::gam( formula(mess), data=dfsub )
        
        ran.x   = range(residuals(gam.x))
        ran.y   = range(residuals(gam.y))
        log_level( INFO, "value=%g.  x residual range = [%g,%g]", value, ran.x[1], ran.x[2] ) 
        log_level( INFO, "value=%g.  y residual range = [%g,%g]", value, ran.y[1], ran.y[2] ) 
        flush.console(); 
        
        #   do not extrapolate if a value is already present
        mask_extrap = Mask3D[ , iV, ]  &  is.na( Lookup4D[ ,iV, ,'x'] )    # mask_extrap is 2D
                      
        for( iC in 2:length(Chroma.vector) )
            {
            chroma  = Chroma.vector[iC]

            mask    = mask_extrap[ ,iC]         # mask is 1D
            
            if( ! any(mask) )   next    # no interpolation needed here
            
            newdata = data.frame( H=Hue.vector[mask], C=chroma )    #V=value,
            
            #   add rectangular coordinates A,B
            tmp = ABfromHC( newdata$H, newdata$C )              
            newdata   = cbind( newdata, A = tmp$A, B = tmp$B )       
                           
            out[mask,iV,iC,'x']  =  predict( gam.x, newdata=newdata )
            out[mask,iV,iC,'y']  =  predict( gam.y, newdata=newdata )
            
            samples = samples + sum(mask)
            }        
        }

    log_level( INFO, "Extrapolated %d samples in %g seconds.", samples, gettime() - time_start )
    
    return( out )
    }
    
    

        
#   returns a logical 3D array
#   for given HVC, the entry is TRUE if an xy value is required for the HVC

makeRequiredMask3D <- function( Munsell2xy )
    {   
    time_start  = gettime()
    
    Hue.vector      = seq(2.5,100,by=2.5)  
            
    Value.vector    = sort( unique(Munsell2xy$V) )
    Value.vector    = c( 0, Value.vector )  # put 0 at beginning, so we can interpolate between 0 and 0.2
    #   Value.vector    = c(0,0.2,0.4,0.6,0.8,1,2,3,4,5,6,7,8,9,10)
    
    #   put 0 at beginning of chroma.vector,
    #   because we will put the neutral x,y there, for Illuminant C
    Chroma.vector   = seq(0,50,by=2) 
    
    nH  = length(Hue.vector)
    nV  = length(Value.vector)
    nC  = length(Chroma.vector)
    
    #   allocate 3D array, and copy Munsell2xy to it
    Mask3D  = array( FALSE, dim=c( nH, nV, nC ) )
    
    dimnames(Mask3D)[[1]] = as.character( Hue.vector )    
    dimnames(Mask3D)[[2]] = as.character( Value.vector )
    dimnames(Mask3D)[[3]] = as.character( Chroma.vector )        
    
    
    if( TRUE )
        #   use the entire all.dat
        dfsub   = Munsell2xy
    else
        #   only the real ones
        dfsub   = Munsell2xy[ Munsell2xy$real,  ]   # ; cat( "nrow(dfsub) = ", nrow(dfsub), '\n' )
        
    
    for( i in 1:nrow(dfsub) )
        {
        row = dfsub[i, ]
        
        iH  = match( row$H, Hue.vector )    # H is always a multiple of 2.5
        iV  = match( row$V, Value.vector )
        iC  = match( row$C, Chroma.vector )  
        
        Mask3D[ iH, iV, iC  ] = TRUE
        }    

        
    #   now dilate with mask
    out = Mask3D
    
    huedelta    = 1     # 2
    chromadelta = 1     # 2
    valuedelta  = 2     # 1 fails when hcinterp='bicub',  vinterp='lin'
    for( iV in 2:nV  )
        {
        mask00  = Mask3D[ , iV, ]
        
        # rotate and translate mask00
        #   hue
        for( h in -huedelta:huedelta )    
            {
            hrot    =  (((1:nH) + h-1) %% nH)  +  1
            maskrot = mask00[ hrot, , drop=FALSE ]
            
            # chroma
            for( j in 0:chromadelta ) 
                {
                cshift  = pmax( (1-j) : (nC-j), 1 )
                mask    = maskrot[  , cshift, drop=FALSE ]
                
                #   copy to valuedelta levels above and below
                vlevel  = max( iV-valuedelta, 2 )  :  min( iV+valuedelta,nV )
                for( v in vlevel )
                    {
                    out[  , v,  ] = mask  |  out[  , v,  ]
                    }
                }
            }
        }

    log_level( INFO, "Made extrapolation mask in %g seconds.",
                            gettime() - time_start )
                            
    count   = integer( nV )                     
    for( iV in 1:nV )
        count[iV]   = sum( out[ ,iV, ] )
        
    count           = c( count, sum(count) )
    names(count)    = c(as.character(Value.vector),'Total')
    print( count )
    
    return( invisible(out) )
    }
    
    