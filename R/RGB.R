


################      Munsell  <-->  RGB    #######################

RGBtoMunsell <- function( RGB, space='sRGB', maxValue=255, adaption='bradford', ... )
    {
    requireNamespace( 'spacesRGB', quietly=TRUE )
    
    theSpace = spacesRGB::getRGB(space)
    if( is.null(theSpace) ) return(NULL)
    
    # Convert RGB input into CIE XYZ and xyY coordinates
    XYZ <- spacesRGB::XYZfromRGB( RGB, space=space, maxValue=maxValue )$XYZ  # for white, Y=100
    if( is.null(XYZ) )  return(NULL)
        
    XYZ = 100 * XYZ
    
    xyY <- XYZ2xyY(XYZ)

    # adapt xyY from the RGB space white to C
    white   = c( theSpace$primaries[4, ], 100 ) 
    white.C = c( p.xyC['NBS',], 100 )

    xyY.adapted = adaptxyY( xyY, white, white.C, method=adaption ) #; print( xyY.adapted )

    # Convert adapted xyY coordinates to Munsell coordinates
    tmp = xyYtoMunsell( xyY.adapted, ... )
    HVC = tmp$HVC
    rownames(HVC) = tmp$SAMPLE_NAME
    HVC
    }

MunsellToRGB <- function( MunsellSpec, space='sRGB', maxValue=255, adaption='bradford', ... )
    {
    requireNamespace( 'spacesRGB', quietly=TRUE )
    
    theSpace = spacesRGB::getRGB(space)
    if( is.null(theSpace) ) return(NULL)
        
    # Convert from Munsell notation to CIE coordiantes
    tmp <- MunsellToxyY( MunsellSpec, ... )
    if( is.null(tmp) )  return(NULL)

    out = tmp
    out$HVC = NULL  # erase this column

    xyY <- tmp$xyY

    # adapt xyY from C to the RGB space white

    white.C = c( p.xyC['NBS',], 100 )    
    white   = c( theSpace$primaries[4, ], 100 )   

    xyY.adapted = adaptxyY( xyY, white.C, white, method=adaption )

    XYZ.adapted <- xyY2XYZ( xyY.adapted )

    # Convert CIE XYZ coordinates to RGB coordinates
    out = cbind( out, spacesRGB::RGBfromXYZ( XYZ.adapted/100, space=space, maxValue=maxValue ) )

    rnames  = tmp$SAMPLE_NAME
    if( anyDuplicated(rnames) ) rnames = 1:nrow(out)

    rownames(out) = rnames
    out
    }




