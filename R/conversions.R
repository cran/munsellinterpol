


################      Munsell  <-->  XYZ    #######################

#   Convert a Munsell specification into XYZ coordinates, by interpolating over the extrapolated Munsell renotation data.
#   returns XYZ that is adapted to C
MunsellToXYZ <- function( MunsellSpec, ... )
    {
    if( ! requireNamespace('spacesXYZ') )   return(NULL)
    
    tmp <- MunsellToxyY( MunsellSpec, ... )
    if( is.null(tmp) )  return(NULL)
    
    XYZ <- spacesXYZ::XYZfromxyY( tmp$xyY )
    rownames(XYZ)   = tmp$SAMPLE_NAME
    XYZ
    }
    
#   convert XYZ to an HVC matrix
#   XYZ     must be already adapted to C
XYZtoMunsell <- function( XYZ, ... )
    {
    if( ! requireNamespace('spacesXYZ') )   return(NULL)
    
    xyY = spacesXYZ::xyYfromXYZ(XYZ)
    if( is.null(xyY) )  return(NULL)
    
    tmp = xyYtoMunsell( xyY, ... )
    HVC = tmp$HVC
    
    rnames  = rownames(XYZ)
    if( is.null(rnames) )   rnames = tmp$SAMPLE_NAME
    rownames(HVC)   = rnames
    
    HVC
    }




################      Munsell  <-->  Lab    #######################
#   Convert a Munsell specification into CIE Lab coordinates, by interpolating over the extrapolated Munsell renotation data.
MunsellToLab <- function( MunsellSpec, white=c(95.047,100,108.883), adapt='Bradford', ... )
    {
    if( ! requireNamespace('spacesXYZ') )   return(NULL)
    
    XYZ = MunsellToXYZ( MunsellSpec, ... )
    if( is.null(XYZ) )  return(NULL)

    #   adapt this XYZ from C to given white
    white.C     = spacesXYZ::XYZfromxyY( c( p.xyC['NBS',],100) )

    theCAT      = spacesXYZ::CAT( white.C, white, method=adapt )
        
    XYZ.adapted = spacesXYZ::adaptXYZ( theCAT, XYZ )

    Lab         = spacesXYZ::LabfromXYZ( XYZ.adapted, white )

    rownames(Lab)   = rownames(XYZ)
    
    Lab
    }

LabtoMunsell <- function( Lab, white=c(95.047,100,108.883), adapt='Bradford', ... )
    {
    if( ! requireNamespace('spacesXYZ') )   return(NULL)

    # make CAT from given white to illuminant C    
    white.C = spacesXYZ::XYZfromxyY( c( p.xyC['NBS',],100) )   
    
    theCAT  = spacesXYZ::CAT( white, white.C, method=adapt )
    if( is.null(theCAT) )   return(NULL)
    
    XYZ = spacesXYZ::XYZfromLab( Lab, white )

    XYZ.adapted = spacesXYZ::adaptXYZ( theCAT, XYZ )
    
    if( is.null(XYZ.adapted) )  return(NULL)

    HVC = XYZtoMunsell( XYZ.adapted, ... )
    #HVC = tmp$HVC
    #rownames(HVC)   = tmp$SAMPLE_NAME
    return(HVC)
    }




################      Munsell  <-->  Luv    #######################
#Convert a Munsell specification into CIE Luv coordinates, by interpolating over the extrapolated Munsell renotation data.

MunsellToLuv <- function( MunsellSpec, white=c(95.047,100,108.883), adapt='Bradford', ... )
    {
    if( ! requireNamespace('spacesXYZ') )   return(NULL)
    
    tmp = MunsellToxyY( MunsellSpec, ... )
    if( is.null(tmp) )  return(NULL)

    XYZ = spacesXYZ::XYZfromxyY( tmp$xyY )

    # make CAT from illuminant C  to   given white 
    white.C = spacesXYZ::XYZfromxyY( c( p.xyC['NBS',],100) )   
    
    theCAT  = spacesXYZ::CAT( white.C, white, method=adapt )
    if( is.null(theCAT) )   return(NULL)
        
    # adapt XYZ from C to given white
    XYZ.adapted = spacesXYZ::adaptXYZ( theCAT, XYZ )

    Luv = spacesXYZ::LuvfromXYZ( XYZ.adapted, white )
    rownames(Luv)   = tmp$SAMPLE_NAME
    
    return(Luv)
    }

LuvtoMunsell <- function( Luv, white=c(95.047,100,108.883), adapt='Bradford', ... )
    {
    if( ! requireNamespace('spacesXYZ') )   return(NULL)
        
    XYZ = spacesXYZ::XYZfromLuv( Luv, white )
    
    # make CAT from given white to illuminant C    
    white.C = spacesXYZ::XYZfromxyY( c( p.xyC['NBS',],100) )   
    
    theCAT  = spacesXYZ::CAT( white, white.C, method=adapt )
    if( is.null(theCAT) )   return(NULL)    

    XYZ.adapted = spacesXYZ::adaptXYZ( theCAT, XYZ  )

    HVC <- XYZtoMunsell( XYZ.adapted, ... )
    HVC
    }




################      Munsell  <-->  sRGB    #######################
# Convert sRGB input into CIE XYZ and xyY coordinates

sRGBtoMunsell <- function( sRGB, maxSignal=255, ... )
    {
    if( ! requireNamespace('spacesXYZ') )   return(NULL)

    if( ! requireNamespace('spacesRGB') )   return(NULL)

    XYZ = spacesRGB::XYZfromRGB( sRGB, space='sRGB', which='scene', maxSignal=maxSignal )$XYZ 
    if( is.null(XYZ) )  return(NULL)

    # adapt XYZ from D65 to C, using precomputed p.D65toC_CAT    
    XYZ.adapted = spacesXYZ::adaptXYZ( p.D65toC_CAT, XYZ ) #; print( xyY.adapted )

    # Convert adapted XYZ coordinates to Munsell coordinates
    HVC = XYZtoMunsell( XYZ.adapted, ... )
    if( is.null(HVC) )  return(NULL)

    # rownames(HVC) = tmp$SAMPLE_NAME
    
    HVC
    }


MunsellTosRGB <- function( MunsellSpec, maxSignal=255, ... )
    {
    if( ! requireNamespace('spacesXYZ') )   return(NULL)

    if( ! requireNamespace('spacesRGB') )   return(NULL)
    
    # Convert from Munsell notation to CIE coordiantes
    tmp = MunsellToxyY( MunsellSpec, ... )
    if( is.null(tmp) )  return(NULL)

    out = tmp
    out$HVC = NULL  # erase this column

    xyY = tmp$xyY

    # adapt xyY from C to D65, using precomputed p.CtoD65_CAT 
    xyY.adapted = spacesXYZ::adaptxyY( p.CtoD65_CAT, xyY )
    
    XYZ.adapted = spacesXYZ::XYZfromxyY( xyY.adapted )

    # Convert CIE XYZ coordinates to sRGB coordinates
    out = cbind( out, spacesRGB::RGBfromXYZ( XYZ.adapted, space='sRGB', which='scene', maxSignal=maxSignal ) )

    rnames  = tmp$SAMPLE_NAME
    if( anyDuplicated(rnames) ) rnames = 1:nrow(out)

    rownames(out) = rnames
    out
    }


#   xyz2srgb is here only because package 'colorscience' depends on it
xyz2srgb <- function(XYZ){
if ((length(XYZ) %% 3) != 0)  stop('XYZ matrix must be n x 3')
if (is.null(dim(XYZ))) if (length(XYZ)>2) XYZ<-matrix(XYZ, ncol=3,byrow=TRUE)
M <- matrix(c(3.2406, -1.5372, -0.4986, -0.9689, 1.8758, 0.0415, 0.0557, -0.2040, 1.0570),3,3,byrow=TRUE)
RGB <- t(M %*% t(XYZ))
# START: lines added March 2013 to set out-of-gamut flag.  
# The out-of-gamut flag is a column vector of Boolean true/false values.  Each
# entry corresponds to one row of the input matrix XYZ.
NumberOfInputs <- dim(RGB)[1]
OutOfGamutFlag <- -99 * matrix(1,NumberOfInputs)
for (index in 1:NumberOfInputs){
if (RGB[index,1] < 0 || RGB[index,1] > 1 || RGB[index,2] < 0 || RGB[index,2] > 1 || RGB[index,3] < 0 || RGB[index,3] > 1) OutOfGamutFlag[index] = TRUE else OutOfGamutFlag[index] <- FALSE
}
# END: lines added March 2013 to set out-of-gamut flag
RGB[which(RGB<0)] <- 0
RGB[which(RGB>1)] <- 1
DACS = matrix(0,dim(XYZ)[1],3)
index <- RGB<=0.0031308
index[index==TRUE] <- 1
index[index==FALSE] <- 0
DACS <- DACS+index*(12.92*RGB)
DACS <- DACS+(1-index)*(1.055*RGB^(1/2.4)-0.055)
RGB <- ceiling(DACS*255)
RGB[which(RGB<0)] <- 0
RGB[which(RGB>255)] <- 255
return(list(Status.ind  = 1, sRGB=RGB, OutOfGamutFlag=OutOfGamutFlag))
}

srgb2xyz <- function(RGBmatrix){
if ((length(RGBmatrix) %% 3) != 0)  stop('RGBmatrix must be n x 3')
if (is.null(dim(RGBmatrix))) if (length(RGBmatrix)>2) RGBmatrix<-matrix(RGBmatrix, ncol=3,byrow=TRUE)
XYZmatrix <- matrix(0,dim(RGBmatrix)[1],3)
M <- matrix(c(0.4124, 0.3576, 0.1805, 0.2126, 0.7152, 0.0722, 0.0193, 0.1192, 0.9505),3,3,byrow=TRUE)
DACS <- RGBmatrix/255
RGBmatrix <- matrix(0,dim(DACS)[1],3)
index <- DACS <= 0.04045
index[index==TRUE] <- 1
index[index==FALSE] <- 0
RGBmatrix = RGBmatrix + (index)*(DACS/12.92)
RGBmatrix = RGBmatrix + (1-index)*((DACS+0.055)/1.055)^2.4
XYZmatrix <- t(M %*% t(RGBmatrix))
XYZmatrix
}


#############     Lightness in [0,100]  <-->  linear Y in [0,1]   ###############
#
#   these threshold=(24/116)^3 and coeff=(12/116)^3    are taken from CIE 15:2004 Section 8.2
#   see http://brucelindbloom.com/LContinuity.html.  Glenn Davis

#   linear Y    in the interval [0,1], nominally
#   returns L in [0,100]
#
#   Some of the x,y in all.dat are negative, and sometimes x+y > 1, which implies z<-.
#   So for Munsell purposes we must allow the tristimulus coordinates XYZ to be negative.
#   To accomodate this, Lightness_from_linear() is defined for negative numbers
#   so that it is an odd function.

Lightness_from_linear  <-  function( Y )
    {
    #   thresh = (24/116)^3
    
    s = sign(Y)    
    Y = s * Y       # Y is now non-negative
    
    out = ifelse( Y < (24/116)^3, (116/12)^3 * Y, 116*Y^(1/3) - 16 )
    
    #   correct any negatives by multiplying by s
    return( s * out )
    }
    
#   nonlinear L         in the interval [0,100]    
#   returns linear Y    in the interval [0,1]
linear_from_Lightness  <-  function( L )
    {    
    s   = sign(L)
    L   = s * L     # L is now non-negative
    
    out = ifelse( L < 8, (12/116)^3 * L, ((L + 16)/116)^3 )
    
    #   correct any negatives by multiplying by s
    return( s * out )    
    }
    
    

        