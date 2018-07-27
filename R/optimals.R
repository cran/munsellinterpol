


IsWithinMacAdamLimits <- function( xyY, Illuminant='C' )
    {
    if( ! requireNamespace( 'geometry', quietly=TRUE ) )    return( NULL )
    
    if( ! requireNamespace( 'spacesXYZ', quietly=TRUE ) )   return( NULL )
        
        
        
    xyY = prepareNx3( xyY )    
    if( is.null(xyY) )  return(NULL)
    
    full    = names( p.OptimalHull )
    idx     = pmatch( toupper(Illuminant), full )
    if( is.na(idx) )
        {
        log.string( ERROR, "Illuminant='%s' is invalid.", Illuminant )
        return(NULL)
        }
    Illuminant  = full[idx]
    
    hull    = p.OptimalHull[[ Illuminant ]]

    out = is.finite( geometry::tsearchn( hull$XYZ, hull$tessellation, spacesXYZ::XYZfromxyY(xyY) )$idx )
    
    return( out )
    }
    

#   xyY     Nx3 matrix of points on the optimal boundary
#
#   returns a list with 2 components
#       XYZ             the given points xyY -> XYZ, an Nx3 matrix
#       tesselation     as returned from geometry::delaunayn
#
makeOptimalHull <- function( xyY )
    {
    if( ! requireNamespace( 'spacesXYZ', quietly=TRUE ) )   return( NULL )
    
    if( ! requireNamespace( 'geometry', quietly=TRUE ) )    return( NULL )    
    
    XYZ = spacesXYZ::XYZfromxyY( xyY )
    
    #   option QJ is to joggle the input and avoid degenerate simplices
    #   without this option there are 2 degenerate ones for Illuminant C, and these generate a warning message for every query
    tessellation = geometry::delaunayn( XYZ, options='QJ' ) 
    
    out = list( XYZ=XYZ, tessellation=tessellation )
    
    return( out )
    }
    
    
#   Y       for the Y-plane defining the section.  in the interval (0,100)
#   n       number of points in the section to return
#   illum   'C' or 'D65'
#
#   returns a data.frame with n rows and these columns
#       XYZ     an optimal color, which is defined by 
#       lambda0 
#       lambda1
sectionOptimals  <-  function( Y, n=100, illum='C' )
    {
    ok  = (0 < Y)  &  (Y < 100)
    if( ! ok )
        {
        log.string( WARN, "Y=%g is invalid; it must be in interval (0,100).", Y )
        return(NULL)
        }
        
    ok  = illum %in% colnames( p.ACDs )
    if( ! ok )   
        {
        log.string( ERROR, "illum='%s' is invalid.", illum )
        return(NULL)
        }

    
    wave    = p.ACDs[[1]]   # wavelengths are column #1
    
    if( ! identical( wave, p.xyz1931$Wavelength ) )
        {
        log.string( FATAL, "illuminant '%s' and xyz have different wavelengths.", illum )
        return(NULL)
        }
    
    funlist = makeReparamFunctionList( wave, p.ACDs[[illum]], as.matrix(p.xyz1931[ ,2:4]) )
    

    Y0  = ( 0:(n-1) ) / n
    Y1  = (Y0 + Y/100) %% 1
    
    omega0  = funlist$omega.from.integral[[2]]( Y0 )
    omega1  = funlist$omega.from.integral[[2]]( Y1 )
    
    bandstop    = (omega1 < omega0)
    
    XYZ.white   = numeric(3)
    XYZ =   matrix( NA_real_, n, 3 )
    
    for( j in 1:3 )
        {
        XYZ.white[j]    = 100 * funlist$integral.from.omega[[j]]( 1 )
        XYZ[ ,j]        = 100 * ( funlist$integral.from.omega[[j]]( omega1 ) - funlist$integral.from.omega[[j]]( omega0 ) )
        
        #   correct the bandstops
        XYZ[bandstop,j] = XYZ[bandstop,j]  +  XYZ.white[j]
        }
        
    #   force Y to be exact
    XYZ[ ,2]    = Y
        
    out         = data.frame( row.names=1:n )
    out$XYZ     = XYZ
    out$lambda0 = funlist$lambda.from.omega( omega0 )
    out$lambda1 = funlist$lambda.from.omega( omega1 )
    
    return(out)
    }
    
    