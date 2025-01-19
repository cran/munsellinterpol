


################      Munsell  <-->  RGB    #######################

RGBtoMunsell <- function( RGB, space='sRGB', which='scene', maxSignal=255, adapt='Bradford', ... )
    {
    for( p in c('spacesRGB','spacesXYZ') )
        {
        if( ! requireNamespace( p, quietly=TRUE ) )
            {
            log_level( ERROR, "required package '%s' could not be loaded.", p )
            return(NULL)
            }
        }
    
    #theSpace = spacesRGB::getRGB(space,full=FALSE)
    #if( is.null(theSpace) ) return(NULL)
    
    RGB = prepareNx3( RGB )
    if( is.null(RGB) )  return(NULL)    
    
    # Convert RGB input into CIE XYZ and xyY coordinates
    XYZ = spacesRGB::XYZfromRGB( RGB, space=space, which=which, maxSignal=maxSignal )$XYZ  # for white, Y=100
    if( is.null(XYZ) )  return(NULL)
    
    #   xyY = spacesXYZ::xyYfromXYZ(XYZ)

    # adapt xyY from the RGB space white to C
    white       = spacesRGB::getWhiteXYZ( space, which=which )  # theSpace$whiteXYZ
    
    white.C     = spacesXYZ::XYZfromxyY( c( p.xyC['NBS',], 100 ) )

    theCAT      = spacesXYZ::CAT( white, white.C, method=adapt )
    
    #   xyY.adapted = spacesXYZ::adaptxyY( theCAT, xyY )   #; print( xyY.adapted )

    XYZ.adapted = spacesXYZ::adaptXYZ( theCAT, XYZ ) #; print( xyY.adapted )
    
    # Convert adapted xyY coordinates to Munsell coordinates
    #tmp = xyYtoMunsell( xyY.adapted, ... )
    #HVC = tmp$HVC
    #rownames(HVC) = tmp$SAMPLE_NAME
    
    # Convert adapted XYZ coordinates to Munsell coordinates
    HVC = XYZtoMunsell( XYZ.adapted, ... )
    if( is.null(HVC) )  return(NULL)

    #   do an additional test for neutrals
    #   to correct possible roundoff error
    neutral = RGB[ ,1]==RGB[ ,2]    &    RGB[ ,2]==RGB[ ,3]

    neutral[ is.na(neutral) ]   = FALSE

    if( any(neutral) )
        {
        # when both R==G==B, then set H=0 and C=0
        HVC[neutral,c(1,3)]     = 0

        if( is.null(rownames(RGB)) )
            #   RGB has no rownames, which means that rownames(HVC) is MunsellName
            rownames(HVC)[neutral]  = MunsellNameFromHVC( HVC[neutral, ] )
        }
    
    return( HVC )
    }

MunsellToRGB <- function( MunsellSpec, space='sRGB', which='scene', maxSignal=255, adapt='Bradford', ... )
    {
    for( p in c('spacesRGB','spacesXYZ') )
        {
        if( ! requireNamespace( p, quietly=TRUE ) )
            {
            log_level( ERROR, "required package '%s' could not be loaded.", p )
            return(NULL)
            }
        }
    
    #theSpace = spacesRGB::getRGB(space,full=FALSE)
    #if( is.null(theSpace) ) return(NULL)
        
    # Convert from Munsell notation to CIE coordiantes
    tmp = MunsellToxyY( MunsellSpec, ... )
    if( is.null(tmp) )  return(NULL)

    out = tmp
    
    #   look for chroma==0 exactly
    neutral = out$HVC[ ,3] == 0
    neutral[ is.na(neutral) ]   = FALSE    
    
    out$HVC = NULL  # erase this column

    #xyY = tmp$xyY
    
    XYZ = spacesXYZ::XYZfromxyY( out$xyY )
    
    # adapt xyY from C to the RGB space white
    white.C =   spacesXYZ::XYZfromxyY( c( p.xyC['NBS',], 100 ) )    
    
    white   =   spacesRGB::getWhiteXYZ( space, which=which )     #theSpace$whiteXYZ
    if( is.null(white) )    return(NULL)
    
    theCAT      = spacesXYZ::CAT( white.C, white, method=adapt )
    
    XYZ.adapted = spacesXYZ::adaptXYZ( theCAT, XYZ )

    # XYZ.adapted = spacesXYZ::XYZfromxyY( xyY.adapted )

    # Convert CIE XYZ coordinates to RGB coordinates
    out = cbind( out, spacesRGB::RGBfromXYZ( XYZ.adapted, space=space, which=which, maxSignal=maxSignal ) )

    #rnames  = tmp$SAMPLE_NAME
    #if( anyDuplicated(rnames) ) rnames = 1:nrow(out)
    #rownames(out) = rnames
    
    if( any(neutral) )
        {
        RGB.mean = rowMeans( out$RGB[neutral, , drop=FALSE] )
        
        #   the vector on the right side is replicated over all columns on the left        
        out$RGB[neutral, ]  =  RGB.mean  # cbind( RGB.mean, RGB.mean, RGB.mean )
        }
    
    return( out )
    }




