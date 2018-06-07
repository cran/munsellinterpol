


################      Munsell  <-->  XYZ    #######################

#   returns XYZ that is adapted to C
MunsellToXYZ <- function( MunsellSpec, ... )
{#Convert a Munsell specification into XYZ coordinates, by interpolating over the extrapolated Munsell renotation data.
if( ! is.character(MunsellSpec) ) {
    MunsellSpec = prepareNx3(MunsellSpec)
    if( is.null(MunsellSpec) )  return(NULL)
    }
tmp <- MunsellToxyY( MunsellSpec, ... )
if( is.null(tmp) )  return(NULL)
XYZ <- xyY2XYZ( tmp$xyY )
rownames(XYZ)   = tmp$SAMPLE_NAME
XYZ
}

#   XYZ     must be already adapted to C
XYZtoMunsell <- function( XYZ, ... )
{ # convert XYZ to an HVC matrix
xyY <- XYZ2xyY(XYZ)
if( is.null(xyY) )  return(NULL)
tmp <- xyYtoMunsell( xyY, ... )
HVC = tmp$HVC
rnames  = rownames(XYZ)
if( is.null(rnames) )   rnames = tmp$SAMPLE_NAME
rownames(HVC)   = rnames
HVC
}




################      Munsell  <-->  Lab    #######################

MunsellToLab <- function( MunsellSpec, white=c(95.047,100,108.883), adaption='bradford', ... )
{#Convert a Munsell specification into CIE Lab coordinates, by interpolating over the extrapolated Munsell renotation data.
if( ! is.character(MunsellSpec) ) {
    MunsellSpec = prepareNx3(MunsellSpec)
    if( is.null(MunsellSpec) )  return(NULL)
    }
tmp <- MunsellToxyY(MunsellSpec, ...)
if( is.null(tmp) )  return(NULL)

XYZ <- xyY2XYZ( tmp$xyY )

#   adapt this XYZ from C to white
white.C  = xyY2XYZ( c( p.xyC['NBS',],100) )

XYZ.adapted = adaptXYZ( XYZ, white.C, white, method=adaption )

Lab <- xyz2lab( XYZ.adapted, white )

rownames(Lab)   = tmp$SAMPLE_NAME
return(Lab)
}

labtoMunsell <- function( Lab, white=c(95.047,100,108.883), adaption='bradford', ... )
{
XYZ <- lab2xyz( Lab, white )
if( is.null(XYZ) )  return(NULL)

# adapt XYZ from given white to illuminant C
white.C  = xyY2XYZ( c( p.xyC['NBS',],100) )

XYZ.adapted = adaptXYZ( XYZ, white, white.C, method=adaption )

HVC <- XYZtoMunsell( XYZ.adapted, ... )
#HVC = tmp$HVC
#rownames(HVC)   = tmp$SAMPLE_NAME
return(HVC)
}




################      Munsell  <-->  Luv    #######################

MunsellToLuv <- function( MunsellSpec, white=c(95.047,100,108.883), adaption='bradford', ... )
{#Convert a Munsell specification into CIE Luv coordinates, by interpolating over the extrapolated Munsell renotation data.
tmp <- MunsellToxyY( MunsellSpec, ... )
if( is.null(tmp) )  return(NULL)

XYZ <- xyY2XYZ( tmp$xyY )

# adapt XYZ from C to given white
white.C  = xyY2XYZ( c( p.xyC['NBS',],100) )

XYZ.adapted = adaptXYZ( XYZ, white.C, white, method=adaption )

Luv <- xyz2luv( XYZ.adapted, white )
rownames(Luv)   = tmp$SAMPLE_NAME
return(Luv)
}

luvtoMunsell <- function( Luv, white=c(95.047,100,108.883), adaption='bradford', ... )
{
XYZ <- luv2xyz( Luv, white )

# adapt XYZ from white to illuminant C
# adapt XYZ from given white to illuminant C
white.C  = xyY2XYZ( c( p.xyC['NBS',],100) )

XYZ.adapted = adaptXYZ( XYZ, white, white.C, method=adaption )

HVC <- XYZtoMunsell( XYZ.adapted, ... )
HVC
}




################      Munsell  <-->  sRGB    #######################

sRGBtoMunsell <- function( sRGB, maxValue=255, adaption='bradford', ... ){
# Convert sRGB input into CIE XYZ and xyY coordinates
XYZ <- srgb2xyz(sRGB,maxValue=maxValue)  # for white, Y=100
if( is.null(XYZ) )  return(NULL)

xyY <- XYZ2xyY(XYZ)

# adapt xyY from D65 to C
white.D65   = c( 0.3127, 0.3290, 100 )    # xy are from the official sRGB standard
white.C     = c( p.xyC['NBS',], 100 )

xyY.adapted = adaptxyY( xyY, white.D65, white.C, method=adaption ) #; print( xyY.adapted )

# Convert adapted xyY coordinates to Munsell coordinates
tmp = xyYtoMunsell( xyY.adapted, ... )
HVC = tmp$HVC
rownames(HVC) = tmp$SAMPLE_NAME
HVC
}


MunsellTosRGB <- function( MunsellSpec, maxValue=255, adaption='bradford', ... ){
# Convert from Munsell notation to CIE coordiantes
tmp <- MunsellToxyY( MunsellSpec, ... )
if( is.null(tmp) )  return(NULL)

out = tmp
out$HVC = NULL  # erase this column

xyY <- tmp$xyY

# adapt xyY from C to D65
white.D65   = c( 0.3127, 0.3290, 100 )    # xy are from the official sRGB standard
white.C     = c( p.xyC['NBS',], 100 )

xyY.adapted = adaptxyY( xyY, white.C, white.D65, method=adaption )

XYZ.adapted <- xyY2XYZ( xyY.adapted )

# Convert CIE XYZ coordinates to sRGB coordinates
out = cbind( out, xyz2srgb( XYZ.adapted, maxValue=maxValue ) )

rnames  = tmp$SAMPLE_NAME
if( anyDuplicated(rnames) ) rnames = 1:nrow(out)

rownames(out) = rnames
out
}




################      XYZ  <-->  xyY    #######################

XYZ2xyY <- function( XYZ ){
XYZ = prepareNx3(XYZ)
if( is.null(XYZ) )  return(NULL)

xyY <- cbind(NA_real_, NA_real_, XYZ[ ,2])
rownames(xyY) = rownames(XYZ)
colnames(xyY) = c('x','y','Y')

denom       = rowSums( XYZ )
w <- which(0<denom  &  0<=XYZ[ ,2])
if (length(w)>0){
xyY[w,1]    = XYZ[w,1] / denom[w]    # (XYZ[w,1]+XYZ[w,2]+XYZ[w,3])
xyY[w,2]    = XYZ[w,2] / denom[w]    # (XYZ[w,1]+XYZ[w,2]+XYZ[w,3])
}
xyY
}

xyY2XYZ <- function( xyY ){
xyY = prepareNx3(xyY)
if( is.null(xyY) )  return(NULL)

XYZ <- cbind( NA_real_, xyY[,3], NA_real_)
rownames(XYZ) = rownames(xyY)
colnames(XYZ) = c('X','Y','Z')

w <- which( xyY[,2] != 0  &  0 <= xyY[,3] )    # was 0 < xyY[ ,2]

if (length(w)>0){
xyY_sub = xyY[w, ,drop=FALSE]
mult    =  xyY_sub[ ,3] / xyY_sub[ ,2]
XYZ[w,1] <- mult * xyY_sub[ ,1]
XYZ[w,3] <- mult * (1-xyY_sub[ ,1]-xyY_sub[ ,2])
}

#   treat Y=0 as a special case - pure black
w <- which( xyY[,3] == 0 )
if( length(w) > 0 )
    XYZ[w,1:3] = 0

XYZ
}



################      XYZ  <-->  sRGB    #######################

#   RGBmatrix       n x 3 matrix
#   maxValue        for the device. 255,1023 are common,  Or 1 for sRGB normalized to 1
#
#   the returned XYZ is for D65 white, with Y=100
srgb2xyz <- function( RGBmatrix, maxValue=255 ){
RGBmatrix = prepareNx3(RGBmatrix)
if( is.null(RGBmatrix) )  return(NULL)
                     
RGBmatrix <- RGBmatrix / maxValue

RGBlin = LinearFromDisplay_sRGB( RGBmatrix )

#   matrix with only 4-digit precision, taken from the sRGB standard, and http://en.wikipedia.org/wiki/SRGB  
# M <- matrix(c(0.4124, 0.3576, 0.1805, 0.2126, 0.7152, 0.0722, 0.0193, 0.1192, 0.9505),3,3,byrow=TRUE)

#   matrix with high precision, taken from   http://brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
#M   = matrix( c( 0.4124564,  0.3575761,  0.1804375, 
#                 0.2126729,  0.7151522,  0.0721750, 
#                 0.0193339,  0.1191920,  0.9503041 ), 3, 3, byrow=TRUE )

XYZmatrix <- 100 * tcrossprod( RGBlin, p.sRGB2XYZ )

rownames(XYZmatrix) = rownames(RGBmatrix)
colnames(XYZmatrix) = c('X','Y','Z')
XYZmatrix
}


#   XYZ         adapted to D65 white, with Y=100
#   maxValue    for the device. 1 for RGB normalized to 1
#
#   returns a data.frame with 2 columns
#           RGB         for the device with given maxValue.  clamped to apropriate cube
#           OutOfGamut  logical
xyz2srgb <- function( XYZ, maxValue=255 ){
XYZ = prepareNx3(XYZ)
if( is.null(XYZ) )  return(NULL)

#   matrix with only 4-digit precision, taken from the sRGB standard, and http://en.wikipedia.org/wiki/SRGB    
#M <- matrix(c(3.2406, -1.5372, -0.4986, -0.9689, 1.8758, 0.0415, 0.0557, -0.2040, 1.0570),3,3,byrow=TRUE)

#   matrix with high precision, taken from   http://brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html    
#M   = matrix( c( 3.2404542, -1.5371385, -0.4985314,
#                -0.9692660,  1.8760108,  0.0415560,
#                 0.0556434, -0.2040259,  1.0572252  ), 3, 3, byrow=TRUE )

RGB <- tcrossprod( XYZ/100, p.XYZ2sRGB )     # t(M %*% t(XYZ))    # print( RGB[1, ] - 1 )

# START: lines added March 2013 to set out-of-gamut flag.  
# The out-of-gamut flag is a column vector of Boolean true/false values.  Each
# entry corresponds to one row of the input matrix XYZ.
# There is numerical tolerance here, designed for points in XYZ that should map to vertices of the RGB cube
lo  = -1.e-7
hi  = 1 + 1.e-7
OutOfGamut  = (RGB[ ,1] < lo  |  RGB[ ,1] > hi  |  RGB[ ,2] < lo  |  RGB[ ,2] > hi  |  RGB[ ,3] < lo  |  RGB[ ,3] > hi)
# END: lines added March 2013 to set out-of-gamut flag

#   clamp RGB to the unit cube (inside gamut)
RGB[ RGB<0 ] <- 0
RGB[ RGB>1 ] <- 1

RGBdisplay  = DisplayFromLinear_sRGB( RGB )     # convert from linear to display

RGBdisplay = maxValue * RGBdisplay

colnames(RGBdisplay) = c('R','G','B')

rnames  = rownames(XYZ)
if( is.null(rnames)  ||  0<anyDuplicated(rnames)  )
    #   rnames is no good !  Use trivial names instead.
    rnames = 1:nrow(XYZ)

out = data.frame( row.names=rnames )
out$sRGB        = RGBdisplay
out$OutOfGamut  = OutOfGamut 

return( out )
}




################      XYZ  <-->  Lab    #######################

#   XYZ     Nx3 matrix, or convertible to one
#   white   X,Y,Z of the reference white, for the input XYZ.  Note that the default is D65 with Y=100
xyz2lab <- function( XYZ, white=c(95.047,100,108.883) )
{
XYZ = prepareNx3(XYZ)
if( is.null(XYZ) )  return(NULL)

L   = Lightness_from_linear( XYZ[ ,2]/white[2] )
a   = (500/116) * (Lightness_from_linear( XYZ[ ,1]/white[1] ) - L)
b   = (200/116) * (L - Lightness_from_linear( XYZ[ ,3]/white[3] ))

Lab = cbind(L, a, b)

rownames(Lab) = rownames(XYZ)
colnames(Lab) = c('L','a','b')
Lab
}

#   XYZ     Nx3 matrix, or convertible to one
#   white   X,Y,Z of the reference white, for the output XYZ.  Note that the default is D65 with Y=100.   Observer= 2deg
lab2xyz <- function( Lab, white=c(95.047,100,108.883) )
{
Lab = prepareNx3(Lab)
if( is.null(Lab) )  return(NULL)

X   = white[1] * linear_from_Lightness( (116/500) * Lab[ ,2]  +  Lab[ ,1])
Y   = white[2] * linear_from_Lightness( Lab[ ,1] ) 
Z   = white[3] * linear_from_Lightness(-(116/200) * Lab[ ,3]  +  Lab[ ,1])

XYZ = cbind(X,Y,Z)
rownames(XYZ) = rownames(Lab)
colnames(XYZ) = c('X','Y','Z')

#   treat L==0 as a special case, pure black
w   = which( Lab[ ,1] == 0 )
if( 0 < length(w) )
    XYZ[w, ] = 0

XYZ
}



################      XYZ  <-->  Luv    #######################


#   XYZ     Nx3 matrix, or convertible to one
#   white   X,Y,Z of the reference white, for the input XYZ.  Note that the default is D65 with Y=100 Observer= 2deg
xyz2luv <- function( XYZ, white=c(95.047,100,108.883) ){
#   Charles Poynton. Digital Video and HD. p. 281.
XYZ = prepareNx3(XYZ)
if( is.null(XYZ) )  return(NULL)

#   chromaticity of white
denom   = sum( white * c(1,15,3) )      # white[1] +  15 * white[2]  +  3 * white[3]
ref.u = 4 * white[1] / denom
ref.v = 9 * white[2] / denom

#   chromaticity of XYZs  (1976 UCS)
denom   =  as.numeric( XYZ %*% c(1,15,3) )      # XYZ[,1] + 15 * XYZ[,2]  + 3 * XYZ[,3]
u = 4 * XYZ[,1] / denom
v = 9 * XYZ[,2] / denom

#   fix invalids
bad = (denom <= 0)
if( any(bad) )
    {
    #   these NAs will propagate
    u[bad]  = NA_real_
    v[bad]  = NA_real_
    }

L = Lightness_from_linear( XYZ[ ,2]/white[2] )
u = 13 * L * ( u - ref.u )
v = 13 * L * ( v - ref.v )

Luv = cbind(L,u,v)
rownames(Luv) = rownames(XYZ)
colnames(Luv) = c('L','u','v')

return(Luv)
}

#   Luv     Nx3 matrix, or convertible to one
#   white   X,Y,Z of the reference white, for the output XYZ.  Note that the default is D65 with Y=100 Observer= 2deg
luv2xyz <- function( Luv, white=c(95.047,100,108.883) ){
#   Charles Poynton. Digital Video and HD. p. 281.
Luv = prepareNx3(Luv)
if( is.null(Luv) )  return(NULL)

#   chromaticity of white
denom   = sum( white * c(1,15,3) )      # white[1] +  15 * white[2]  +  3 * white[3]
ref.u = 4 * white[1] / denom
ref.v = 9 * white[2] / denom

u   = Luv[,2] / ( 13 * Luv[ ,1] )  +  ref.u
v   = Luv[,3] / ( 13 * Luv[ ,1] )  +  ref.v

Y   = white[2] * linear_from_Lightness( Luv[ ,1] ) 
X   = Y * (9*u)/(4*v)
Z   = Y * (12 - 3*u - 20*v)/(4*v)

XYZ = cbind(X,Y,Z)
rownames(XYZ) = rownames(Luv)
colnames(XYZ) = c('X','Y','Z')

#   treat L==0 as a special case, pure black
w   = which( Luv[ ,1] == 0 )
if( 0 < length(w) )
    XYZ[w, ] = 0

return(XYZ)
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
    
    
#############     sRGB  EOCF and OECF [0,1]  <-->  [0,1]   ###############    
#   maps [0,1] to [0,1].   Also extended to negative numbers in a symmetric way
#   it is OK if .v is a matrix, and then the return value is a matrix of the same shape
DisplayFromLinear_sRGB <- function( .v )
    {
    s = sign(.v)    
    out = s * .v    #   now non-negative
    out = ifelse( out <= 0.0031308,    12.92 * out,  (1055 * out^(1/2.4) - 55 ) / 1000 )
    return( s * out )
    }

#   maps [0,1] to [0,1].   Also extended to negative numbers in a symmetric way
#   it is OK if .v is a matrix, and then the return value is a matrix of the same shape
LinearFromDisplay_sRGB <- function( .v )
    {
    s = sign(.v)    
    out = s * .v    #   now non-negative
    out = ifelse( out <= 12.92*0.0031308,   out / 12.92,   ( (1000*out + 55)/(1055) ) ^ 2.4  )
    return( s * out )
    }    
    

        