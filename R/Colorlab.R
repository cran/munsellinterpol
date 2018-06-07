

MunsellSpecToColorlabFormat  <-  function( MunsellSpec )
    {
    if( is.character(MunsellSpec) )
        {
        rnames  = MunsellSpec
        MunsellSpec = HVCfromMunsellName( MunsellSpec )
        }
    else
        rnames = NULL
        
    HVC = prepareNx3(MunsellSpec)
    if( is.null(HVC) )  return(NULL)
    
    hue = HVC[ ,1] %% 100
    
    hue.rev = 100 - hue
    
    hue.idx     = as.integer( hue.rev / 10 )
    hue.frac    = hue.rev - 10 * hue.idx
    
    hue.idx     = ((hue.idx-3) %% 10) + 1
    hue.step    = 10 - hue.frac
    
    out = cbind( hue.step, HVC[ ,2:3,drop=FALSE], hue.idx )  #; print( str(out) )
    
    #   fix the neutrals, set both hue number and hue index to 0
    out[ HVC[ ,3]==0, c(1,4) ] = 0
    
    if( is.null(rnames) )   rnames = MunsellNameFromHVC(HVC)
    
    rownames(out)   = rnames
    colnames(out)   = c( 'HS', 'V', 'C', 'HI' )

    return(out)
    }
    
    
ColorlabFormatToMunsellSpec  <-  function( HVCH )    
    {
    HVCH    = prepareNx3( HVCH, M=4 )
    if( is.null(HVCH) ) return(NULL)
    
    hue.rev = 10 * HVCH[ ,4] - HVCH[ ,1]
    
    hue = 100 - ( (hue.rev+30) %% 100 )
    
    out = cbind( hue, HVCH[ ,2:3,drop=FALSE] )
    
    #   fix the neutrals, set hue to 0
    mask.neutral    = (HVCH[ ,3]==0)
    out[  mask.neutral, 1] = 0    
    
    #   invalidate any rows where the hue index is invalid
    mask.ok = (HVCH[ ,4]  %in%  1:10)  |  mask.neutral
    out[ ! mask.ok, ] = NA_real_
    
    rownames(out)   = rownames(HVCH)
    
    if( is.null( rownames(out) ) )
        rownames(out)   = MunsellNameFromHVC(out)
        
    colnames(out)   = c('H','V','C')
    
    return( out )
    }
    