


################      Munsell  <-->  RGB    #######################

RGBtoMunsell <- function( RGB, space='sRGB', maxSignal=255, adapt='Bradford', ... )
    {
    if( ! requireNamespace( 'spacesRGB', quietly=TRUE ) )   return(NULL)
    
    if( ! requireNamespace( 'spacesXYZ', quietly=TRUE ) )   return(NULL)
    
    #theSpace = spacesRGB::getRGB(space,full=FALSE)
    #if( is.null(theSpace) ) return(NULL)
    
    # Convert RGB input into CIE XYZ and xyY coordinates
    XYZ = spacesRGB::XYZfromRGB( RGB, space=space, which='display', maxSignal=maxSignal )$XYZ  # for white, Y=100
    if( is.null(XYZ) )  return(NULL)
        
    #   XYZ = 100 * XYZ
    
    xyY = spacesXYZ::xyYfromXYZ(XYZ)

    # adapt xyY from the RGB space white to C
    white       = spacesRGB::getWhiteXYZ( space, which='display')  # theSpace$whiteXYZ
    white.C     = spacesXYZ::XYZfromxyY( c( p.xyC['NBS',], 100 ) )

    theCAT      = spacesXYZ::CAT( white, white.C, method=adapt )
    
    xyY.adapted = spacesXYZ::adaptxyY( theCAT, xyY )   #; print( xyY.adapted )

    # Convert adapted xyY coordinates to Munsell coordinates
    tmp = xyYtoMunsell( xyY.adapted, ... )
    HVC = tmp$HVC
    rownames(HVC) = tmp$SAMPLE_NAME
    
    HVC
    }

MunsellToRGB <- function( MunsellSpec, space='sRGB', maxSignal=255, adapt='Bradford', ... )
    {
    if( ! requireNamespace( 'spacesRGB', quietly=TRUE ) )   return(NULL)
    
    if( ! requireNamespace( 'spacesXYZ', quietly=TRUE ) )   return(NULL)
    
    #theSpace = spacesRGB::getRGB(space,full=FALSE)
    #if( is.null(theSpace) ) return(NULL)
        
    # Convert from Munsell notation to CIE coordiantes
    tmp = MunsellToxyY( MunsellSpec, ... )
    if( is.null(tmp) )  return(NULL)

    out = tmp
    out$HVC = NULL  # erase this column

    xyY = tmp$xyY

    # adapt xyY from C to the RGB space white
    white.C =   spacesXYZ::XYZfromxyY( c( p.xyC['NBS',], 100 ) )    
    
    white   =   spacesRGB::getWhiteXYZ( space, which='display')     #theSpace$whiteXYZ
    if( is.null(white) )    return(NULL)
    
    theCAT      = spacesXYZ::CAT( white.C, white, method=adapt )
    
    xyY.adapted = spacesXYZ::adaptxyY( theCAT, xyY )

    XYZ.adapted = spacesXYZ::XYZfromxyY( xyY.adapted )

    # Convert CIE XYZ coordinates to RGB coordinates
    out = cbind( out, spacesRGB::RGBfromXYZ( XYZ.adapted, space=space, which='display', maxSignal=maxSignal ) )

    rnames  = tmp$SAMPLE_NAME
    if( anyDuplicated(rnames) ) rnames = 1:nrow(out)

    rownames(out) = rnames
    out
    }




