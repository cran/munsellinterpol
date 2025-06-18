## ----setup, include=FALSE-------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options( width=100 )

## ----echo=TRUE, message=TRUE----------------------------------------------------------------------
library(munsellinterpol)

## ----echo=TRUE, warning=TRUE----------------------------------------------------------------------
load_soil <- function( soil_raw=OxSR::soil_refle, subset=c( "a4", "a6", "a8", "a11", "a14", "a15" ),
                                         wave=380:780 )
    {
    #   locate only the wavelengths in OxSR::soil_refle that are relevant for color
    #   the given data.frame OxSR::soil_refle goes way into the infra-red and these are not needed
    #   We already know that all relevant multiples of 1nm are in soil_raw
    idx     = match( wave, soil_raw$wavelength_nm )

    #   extract just a subset of samples with high variability
    #   the reflectances in OxSR::soil_refle are percentages, so divide by 100    
    mat     = as.matrix(soil_raw)[ idx, subset ] / 100
    
    #   convert the matrix mat to a colorSpec object
    soil_spec = colorSpec::colorSpec( mat, wavelength=wave, quantity='reflectance' )
    
    #   add neutral gray as a test sample
    gray      = colorSpec::neutralMaterial( 0.18, wavelength=wave )
    soil_spec = colorSpec::bind( soil_spec, gray )
   
    return( soil_spec ) 
    }

## ----echo=TRUE, warning=TRUE,  fig.width=7.3, fig.height=5.2,  fig.show='hold'--------------------
soil_spec = load_soil()
par( omi=c(0,0,0,0), mai=c(0.6,0.6,0.3,0.3) )
plot( soil_spec, legend="topleft" )

## ----echo=TRUE, warning=TRUE----------------------------------------------------------------------
soil_data <- function( soil_spec, CCT )
    {
    wave    = colorSpec::wavelength(soil_spec) 
    
    #   make illuminant from the CCT. This is the reference illuminant in IES Standard TM-30.
    illum   = colorSpec::referenceSpectraTM30( CCT, wavelength=wave )
    
    #   make scanner from this illuminant and the standard x,y,z responsivities of the human eye
    scanner = colorSpec::product( illum, "soil", colorSpec::xyz1931.1nm, wavelength=wave )
    
    #   scale so the Perfect Reflecting Diffuser has Y=100; which is conventional for Munsell work
    scanner = colorSpec::calibrate( scanner, response=c(NA,100,NA), method='scaling' )
    
    #   compute white point of illum
    white   = colorSpec::product( colorSpec::neutralMaterial(1,wavelength=wave), scanner )
    
    #   compute XYZ for all soil samples; this is XYZ under illum
    XYZ     = colorSpec::product( soil_spec, scanner )
    
    #   compute Lab for all soil samples, under illum
    Lab = spacesXYZ::LabfromXYZ( XYZ, white )
    
    # chromatically adapt soil sample XYZ under illum, to sample XYZ under Illuminant C
    white.C = 100 * spacesXYZ::standardXYZ( "C.NBS" )   # use the original NBS variant of this white point
    theCAT  = spacesXYZ::CAT( white, white.C )
    XYZ.C   = spacesXYZ::adaptXYZ( theCAT, XYZ )
    
    #   convert XYZ.C (XYZ under Illuminant C) to Munsell HVC.
    HVC     = XYZtoMunsell( XYZ.C, xyC="NBS" )  #  use the original NBS variant, as in standardXYZ() above
    rownames(HVC)   = colorSpec::specnames(soil_spec)
    
    # create output data.frame and add 4 columns: HVC, Munsell notation, ISCC-NBS color name, and Lab
    out = data.frame( row.names=rownames(HVC) )
    out$HVC     = HVC
    out$Munsell = MunsellNameFromHVC( HVC )
    out[[ "ISCC-NBS Name" ]] = ColorBlockFromMunsell( HVC )$Name
    out$Lab     = Lab
    
    return( out )
    }

## ----echo=TRUE, message=TRUE----------------------------------------------------------------------
dat_soil = soil_data( soil_spec, CCT=6500 )
dat_soil

## ----echo=TRUE, message=TRUE----------------------------------------------------------------------
dat_rnd     = roundHVC( dat_soil$HVC, books="soil" )
dat_rnd$HVC = NULL     # delete matrix HVC, because it is redundant in this context
dat_rnd

## ----echo=TRUE, message=TRUE, fig.width=7.3, fig.height=5.2,  fig.show='hold'---------------------
par( omi=c(0,0,0,0), mai=c(0.4,0.4,0.25,0) )
plotPatchesH( "2.5YR",   value=c(2.5,3:8), chroma=c(1,2,3,4,6,8) )
text( c(4,4,3), c(5,4,4), c("a6","a8","a15"), adj=c(0.5,0.5), col="white" )

## ----echo=TRUE, message=TRUE----------------------------------------------------------------------
tbl = data.frame( row.names=1:(2*nrow(dat_soil)) )
tbl[[ "Sample" ]]   = as.character( rbind(rownames(dat_soil),'') )
tbl[[ "Precise" ]] = as.character( rbind('',MunsellNameFromHVC( dat_soil$HVC, digits=3 )) )
tbl[[ "Rounded" ]]  = as.character( rbind('',dat_rnd$MunsellRounded ) )

library( flextable )
myrt <- regulartable( tbl )
myrt <- width( myrt, j=c(2,3), width=2 )
myrt <- align( myrt, j=1:3, align='center', part='all' )
myrt <- hline( myrt, i=seq( 2, nrow(tbl)-2, by=2 ), border=fp_border_default( color="black", width=2 ) )
myrt <- hrule( myrt, rule = "exact" )
idxcolor <- seq(1,nrow(tbl)-1,by=2)
myrt <- height( myrt, i=idxcolor, height=1 )
myrt <- bg( myrt, i=idxcolor, j=2, bg=rgb( round(MunsellTosRGB(dat_soil$HVC)$RGB), max=255 ) )
myrt <- bg( myrt, i=idxcolor, j=3, bg=rgb( round(MunsellTosRGB(dat_rnd$MunsellRounded)$RGB), max=255 ) )
myrt

## ----echo=FALSE, results='asis'-------------------------------------------------------------------
sessionInfo()

