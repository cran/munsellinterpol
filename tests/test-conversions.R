
#   This is for all EXCEPT   HVC  <--->  xyY   transforms

library( munsellinterpol )

options( width=144 )
    
printf <- function( msg, ... )
    {    
    mess = sprintf( msg[1], ... )    # should this really be msg[1] ?
    cat( mess, '\n' )   #, file=stderr() )
    }
    
#   returns time in seconds, from an arbitrary origin
gettime <- function()
    {
    if( requireNamespace('microbenchmark') )
        return( microbenchmark::get_nanotime() * 1.e-9 )
    else
        return( as.double( base::Sys.time() ) )
    }
    
    
testXYZ <- function()
    {
    printf( "---------------------  testXYZ()  -----------------------" )
    require( spacesRGB )
        
    # make random XYZs    
    set.seed(0)    
    count   = 1000
    
    XYZ = matrix( abs(rnorm(3*count)), ncol=3 )
    rownames(XYZ)   = sprintf( "%04d", 1:count )
    
    #---  xyY ---##
    time_start      = gettime()    
    xyY             = XYZ2xyY( XYZ )
    XYZ.back        = xyY2XYZ( xyY )
    time_elapsed    = gettime() - time_start
    
    delta   = rowSums( abs(XYZ - XYZ.back) ) 
    printf( "XYZ -> xyY -> XYZ    max(delta)=%g   %d samples at %g sec/sample",
                            max(delta), count, time_elapsed/count )
    
    failures = sum( 5.e-15 < delta )
    
    if( 0 < failures )
        {
        idx = which.max(delta)        
        printf( "There were %d  XYZ -> xyY -> XYZ failures.  Max error = %g",
                    failures, delta[idx] )

        df  = data.frame( row.names=1 )
        df$XYZ      = XYZ[idx, ,drop=FALSE]
        df$xyY      = xyY[idx, ,drop=FALSE]
        df$XYZ.back = XYZ.back[idx, ,drop=FALSE]        
        print( df )
        
        return(FALSE)
        }
        
    #   test pure black
    black   = c(0,0,0)
    if( ! identical( black, as.numeric(xyY2XYZ( XYZ2xyY(black) ) ) ) )
        {
        printf(  "XYZ -> xyY -> XYZ.back .  pure black not preserved." )
        return(FALSE)
        }
        
    #   test rownames
    if( ! identical( rownames(XYZ), rownames(XYZ.back) ) )
        {
        printf(  "XYZ -> xyY -> XYZ.back .  rownames not preserved." )
        return(FALSE)
        }

        
        
    #---  Lab ---##
    time_start      = gettime()        
    Lab         = xyz2lab( XYZ )
    XYZ.back    = lab2xyz( Lab )
    time_elapsed    = gettime() - time_start
    
    delta   = rowSums( abs(XYZ - XYZ.back) )        
    
    printf( "XYZ -> Lab -> XYZ    max(delta)=%g   %d samples at %g sec/sample", 
                    max(delta), count, time_elapsed/count )
        
    failures = sum( 1.e-14 < delta )
    
    if( 0 < failures )
        {
        idx = which.max(delta)        
        printf(  "There were %d  XYZ -> Lab -> XYZ failures.  Max error = %g",
                    failures, delta[idx] )

        df  = data.frame( row.names=1 )
        df$XYZ      = XYZ[idx, ,drop=FALSE]
        df$Lab      = Lab[idx, ,drop=FALSE]
        df$XYZ.back = XYZ.back[idx, ,drop=FALSE]        
        print( df )
        
        return(FALSE)
        }
        
    #   test pure black
    black   = c(0,0,0)
    if( ! identical( black, as.numeric(lab2xyz( xyz2lab(black) ) ) ) )
        {
        printf(  "XYZ -> Lab -> XYZ.back .  pure black not preserved." )
        return(FALSE)
        }
                        
    #   test rownames
    if( ! identical( rownames(XYZ), rownames(XYZ.back) ) )
        {
        printf(  "XYZ -> Lab -> XYZ.back .  rownames not preserved." )
        return(FALSE)
        }

        
    #---  Luv ---##
    time_start      = gettime()           
    Luv             = xyz2luv( XYZ )
    XYZ.back        = luv2xyz( Luv )
    time_elapsed    = gettime() - time_start
    
    delta   = rowSums( abs(XYZ - XYZ.back) )      
    
    printf( "XYZ -> Luv -> XYZ    max(delta)=%g   %d samples at %g sec/sample", 
                        max(delta), count, time_elapsed/count )
            
    failures = sum( 1.e-13 < delta )
    
    if( 0 < failures )
        {
        idx = which.max(delta)        
        printf(  "There were %d  XYZ -> Luv -> XYZ failures.  Max error = %g",
                    failures, delta[idx] )

        df  = data.frame( row.names=1 )
        df$XYZ      = XYZ[idx, ,drop=FALSE]
        df$Luv      = Luv[idx, ,drop=FALSE]
        df$XYZ.back = XYZ.back[idx, ,drop=FALSE]        
        print( df )
        
        return(FALSE)
        }
        
    #   test pure black
    black   = c(0,0,0)
    if( ! identical( black, as.numeric(luv2xyz( xyz2luv(black) ) ) ) )
        {
        printf( "XYZ -> Luv -> XYZ.back .  pure black not preserved." )
        return(FALSE)
        }
        
    #   test rownames
    if( ! identical( rownames(XYZ), rownames(XYZ.back) ) )
        {
        printf(  "XYZ -> Luv -> XYZ.back .  rownames not preserved." )
        return(FALSE)
        }
        

        
        
    #---  sRGB ---##
    RGB         = matrix( runif(3*count,max=255), ncol=3 )
    rownames(RGB)   = sprintf( "%04d", 1:count )
    
    time_start      = gettime()               
    XYZ             = srgb2xyz( RGB )    
    RGB.back        = xyz2srgb( XYZ )$sRGB      #; print( 'RGB OK' )
    time_elapsed    = gettime() - time_start
    
    delta   = rowSums( abs(RGB - RGB.back) )  
    
    printf( "sRGB -> XYZ -> sRGB    max(delta)=%g   %d samples at %g sec/sample", 
                        max(delta), count, time_elapsed/count )
     
    tol = 5.e-12   
    failures = sum( tol < delta )   
    if( 0 < failures )
        {
        idx = which.max(delta)        
        printf(  "There were %d  sRGB -> XYZ -> sRGB failures.  Max error = %g",
                    failures, delta[idx] )

        df  = data.frame( row.names=1 )
        df$sRGB         = RGB[idx, ,drop=FALSE]
        df$XYZ          = XYZ[idx, ,drop=FALSE]
        df$sRGB.back    = RGB.back[idx, ,drop=FALSE]        
        print( df )
        
        return(FALSE)
        }
        
    #   test pure black
    black   = c(0,0,0)
    if( ! identical( black, as.numeric(xyz2srgb( srgb2xyz(black) )$sRGB ) ) )
        {
        printf(  "sRGB -> XYZ -> sRGB.back .  pure black not preserved." )
        return(FALSE)
        }        
    #   test rownames
    if( ! identical( rownames(RGB), rownames(RGB.back) ) )
        {
        printf( "RGB -> xyY -> RGB.back .  rownames not preserved." )
        return(FALSE)
        }

        
    #---  AdobeRGB ---##
    time_start      = gettime()        
    XYZ         = spacesRGB::XYZfromRGB( RGB, space='AdobeRGB', max=255 )$XYZ
    RGB.back    = spacesRGB::RGBfromXYZ( XYZ, space='AdobeRGB', max=255 )$RGB      #; print( 'RGB OK' )
    time_elapsed    = gettime() - time_start
    
    delta   = rowSums( abs(RGB - RGB.back) )  
    
    printf( "AdobeRGB -> XYZ -> AdobeRGB    max(delta)=%g   %d samples at %g sec/sample", 
                        max(delta), count, time_elapsed/count )
     
    tol = 5.e-10  
    failures = sum( tol < delta )   
    if( 0 < failures )
        {
        idx = which.max(delta)        
        printf(  "There were %d  AdobeRGB -> XYZ -> AdobeRGB failures.  Max error = %g > %g",
                    failures, delta[idx], tol )

        df  = data.frame( row.names=1 )
        df$RGB          = RGB[idx, ,drop=FALSE]
        df$XYZ          = XYZ[idx, ,drop=FALSE]
        df$sRGB.back    = RGB.back[idx, ,drop=FALSE]        
        print( df )
        
        return(FALSE)
        }        
        
    return( TRUE )
    }
     
     
     
testMunsell <- function()
    {
    printf( "---------------------  testMunsell()  -----------------------" )
            
    set.seed(0)    
        
    # make random RGBs in the cube       
    count   = 1000
    RGB = matrix( runif(3*count,min=0,max=255), ncol=3 )
    rownames(RGB)   = sprintf( "%04d", 1:count )
         
    printf( "testing sRGB -> HVC -> sRGB with %d samples...", count )
                
    time_start  = gettime()    
    HVC         = sRGBtoMunsell( RGB )          #; print( str(HVC) )
    RGB.back    = MunsellTosRGB( HVC )$sRGB     #; print( str(RGB.back) )
    time_elapsed    = gettime() - time_start    
    
    mask    = is.na( RGB.back[ ,1] )
    if( any( mask ) )
        {
        printf( "There were %d  sRGB -> HVC -> sRGB  failures (out of %d).", sum(mask), count )
        return(FALSE)
        }
        
    delta   = rowSums( abs(RGB - RGB.back) )  #; print( max(delta) )
    
    printf( "sRGB  max(delta)=%g   (out of %d)",  max(delta,na.rm=TRUE), sum(is.finite(delta)) )
        
    tol     = 2.e-2
    failures = sum( tol < delta ) #; print( failures )
    
    if( 0 < failures )
        {
        idx = which.max(delta)        
        printf( "There were %d  RGB -> HVC -> RGB tolerance failures (out of %d).  Max error = %g > %g.",
                    failures, count, delta[idx], tol )

        df  = data.frame( row.names=1 )
        df$RGB      = RGB[idx, ,drop=FALSE]
        df$HVC      = HVC[idx, ,drop=FALSE]
        df$RGB.back = RGB.back[idx, ,drop=FALSE]        
        print( df )
        
        return(FALSE)
        }
            
    printf( "passed.  %d round trips in %g sec,  %g/sample", count, time_elapsed, time_elapsed/count )
                        
        
    #   test pure black
    black       = c(0,0,0)
    black.back  = as.numeric( MunsellTosRGB( sRGBtoMunsell(black) )$sRGB )
    if( ! identical( black, black.back ) )
        {
        printf( "RGB -> HVC -> RGB .  pure black not preserved." )
        return(FALSE)
        }     

        
    #---   AdobeRGB  -----#
        
    printf( "testing AdobeRGB -> HVC -> AdobeRGB with %d samples...", count )
                
    time_start      = gettime()    
    HVC             = RGBtoMunsell( RGB, space='AdobeRGB' )          #; print( str(HVC) )
    RGB.back        = MunsellToRGB( HVC, space='AdobeRGB' )$RGB     #; print( str(RGB.back) )
    time_elapsed    = gettime() - time_start    
    
    mask    = is.na( RGB.back[ ,1] )
    if( any( mask ) )
        {
        printf( "There were %d  AdobeRGB -> HVC -> AdobeRGB  failures (out of %d).", sum(mask), count )
        return(FALSE)
        }
        
    delta   = rowSums( abs(RGB - RGB.back) )  #; print( max(delta) )
    
    printf( "AdobeRGB  max(delta)=%g   (out of %d)",  max(delta,na.rm=TRUE), sum(is.finite(delta)) )
                
    tol     = 0.5
    failures = sum( tol < delta ) #; print( failures )
    
    if( 0 < failures )
        {
        idx = which.max(delta)        
        printf( "There were %d  AdobeRGB -> HVC -> AdobeRGB tolerance failures (out of %d).  Max error = %g > %g.",
                    failures, count, delta[idx], tol )
        df  = data.frame( row.names=1 )
        df$RGB      = RGB[idx, ,drop=FALSE]
        df$HVC      = HVC[idx, ,drop=FALSE]
        df$RGB.back = RGB.back[idx, ,drop=FALSE]        
        print( df )
        return(FALSE)
        }
            
    printf( "passed.  %d round trips in %g sec,  %g/sample", count, time_elapsed, time_elapsed/count )
                  
        


    #---   Lab   ----#
    printf( "testing Lab -> HVC -> Lab with %d samples...", count )

    time_start      = gettime()        
    Lab             = MunsellToLab( HVC )    
    Lab.back        = MunsellToLab( labtoMunsell( Lab ) )
    time_elapsed    = gettime() - time_start      
    
    if( any( is.na(Lab.back) ) )
        {
        printf( "There were Lab -> HVC  failures (out of %d).", count )        
        return(FALSE)
        }
    
    delta   = rowSums( abs(Lab - Lab.back) )
    
    printf( "Lab  max(delta)=%g   (out of %d)",  max(delta,na.rm=TRUE), sum(is.finite(delta)) )
           
    tol = 1.e-2
    failures = sum( tol < delta, na.rm=TRUE )       
    if( 0 < failures )
        {
        idx = which.max(delta)        
        printf(  "There were %d  Lab -> HVC -> Lab failures (out of %d).  Max error = %g",
                    failures, sum(is.finite(delta)), delta[idx] )

        df  = data.frame( row.names=1 )
        df$Lab      = Lab[idx, ,drop=FALSE]
        df$HVC      = HVC[idx, ,drop=FALSE]
        df$Lab.back = Lab.back[idx, ,drop=FALSE]        
        print( df )
        
        return(FALSE)
        }        

    printf( "passed.  %d round trips in %g sec,  %g/sample", count, time_elapsed, time_elapsed/count )                        
            
    return( TRUE )        
    }
    
    
testColorlab <- function()
    {
    printf( "---------------------  testColorlab()  -----------------------" )
             
    set.seed(0)    
    count   = 1000
    
    # make random HVCs        
    
    H   = runif( count, min=0, max=100 )
    V   = runif( count, min=0, max=10 )
    C   = runif( count, min=1, max=6 )
    
    HVC = cbind( H, V, C )
    
    # add a few some neutrals
    HVC.n   = cbind( 0, runif(20,min=0,max=10), 0 )
    HVC = rbind( HVC, HVC.n )
    HVC = rbind( HVC, c(0,0,0) )    # add pure black
    
    #   integers are the worst case, so make 1/2 the data integers
    idx = 1:(count/2)
    HVC[idx, ] = round( HVC[idx, ] )

    count   = nrow(HVC)    
    printf( "testing HVC -> HVCH -> HVC with %d samples...", count )
        
    time_start    = gettime()
    
    HVCH        = MunsellSpecToColorlabFormat( HVC )
    HVC.back    = ColorlabFormatToMunsellSpec( HVCH )
    
    #   compare HVC and HVC.back 
    HVC.delta   = abs(HVC - HVC.back)
    HVC.delta[ ,1]  = pmin( HVC.delta[ ,1], 100 - HVC.delta[ ,1] )
    delta = rowSums( HVC.delta )    
    
    printf( "Colorlab max(delta) = %g.", max(delta) )
    
    tol = 1.e-16
    failures = sum( tol <= delta )   
    if( 0 < failures )    
        {
        idx = which.max( delta )
        df  = data.frame( row.names=1 )
        df$HVC      = HVC[ idx, , drop=F ]
        df$HVCH     = HVCH[ idx, , drop=F ]
        df$HVC.back = HVC.back[ idx, , drop=F ]
        printf(  "Colorlab maximum inversion error = %g.", delta[idx] )
        print( df )    
        return(FALSE)
        }

    time_elapsed    = gettime() - time_start        
    printf( "passed.  %d round trips in %g sec,  %g/sample", count, time_elapsed, time_elapsed/count )

        
    return( TRUE )
    }
        
x = gettime()   # load microbenchmark
 
if( ! testXYZ() )       stop( "testXYZ() failed !", call.=FALSE )

if( ! testColorlab() )  stop( "testColorlab() failed !", call.=FALSE )

if( ! testMunsell() )   { stop( "testMunsell() failed !", call.=FALSE ) }
   


printf(  "\nPassed all Conversion tests !" )
     