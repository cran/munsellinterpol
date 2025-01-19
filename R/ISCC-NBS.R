#
#   colorBlockFromMunsell()
#
#   MunsellSpec     a numeric Nx3 matrix, with HVC in each row
#                   or a single plain numeric HVC 3-vector
#                   or a character vector of Munsell names
#
#   if MunsellSpec is a 3-vector 
#       MunsellSpec[1]     Munsell hue, on the ASTM D1535 100 point circular scale.  values are wrapped to [1,101)
#       MunsellSpec[2]     Munsell value, must be between 0 and 10
#       MunsellSpec[3]     Munsell chroma, must be non-negative
#
#   if MunsellSpec is an Nx3 matrix
#       each row is defined as above
#
#   if MunsellSpec is an character N-vector of Munsell notations,
#       it is converted to an Nx3 matrix
#
#   return value:
#        data.frame with columns: HVC, ISCC-NBS Number, ISCC-NBS Name
#
#   the function requires these global data.frames:
#       p.System_ISCCNBS    this one is private
#       CentroidsISCCNBS    this one is public
#
#   author:  Glenn Davis


ColorBlockFromMunsell  <-  function( MunsellSpec )
    {    
    if( is.character(MunsellSpec) )
        {
        HVC     = HVCfromMunsellName(MunsellSpec)
        if( is.null(HVC) )  return(NULL)
        
        mask    = is.na( HVC[ ,1] )
        if( any(mask) )
            log_level( ERROR, "%d (out of %d) Munsell notations are invalid.", sum(mask), nrow(HVC) )
        
        rnames  = MunsellSpec
        }
    else  
        {
        HVC = prepareNx3( MunsellSpec )
        if( is.null(HVC) )  return(NULL)

        rnames  = MunsellNameFromHVC( HVC, digits=2 )
        rownames(HVC)   = rnames
        colnames(HVC)   = c('H','V','C')
        }
    
    n   = nrow(HVC)
    if( anyDuplicated(rnames) ) rnames  = 1:n   # rnames will be used later
    
    number      = rep( NA_integer_, n )
    name        = rep( NA_character_, n )
    centroid    = rep( NA_character_, n )
    vmax        = 10                
    
    for( i in 1:n )
        {
        hvc = HVC[i, ]
    
        valid   = all( is.finite(hvc) )  &&  (0 <= hvc[2])  &&  (hvc[2] <= vmax)  &&  (0 <= hvc[3])
        
        if( ! valid )   next
        
        #   wrap hue to the interval [1,101)   (101 is not included)
        hvc[1] = ((hvc[1] - 1 ) %% 100) + 1
            
        if( hvc[2] == vmax )   hvc[2] = vmax - 1.e-6     # because upper comparison below is strict
        
        #   do 6-way comparison.  
        #   Note upper comparisons are strict, and lower comparisons are not strict.
        #   So a point on a boundary is in only 1 block.
        mask.H  = p.System_ISCCNBS$Hmin <= hvc[1]  &  hvc[1] < p.System_ISCCNBS$Hmax
        mask.V  = p.System_ISCCNBS$Vmin <= hvc[2]  &  hvc[2] < p.System_ISCCNBS$Vmax
        mask.C  = p.System_ISCCNBS$Cmin <= hvc[3]  &  hvc[3] < p.System_ISCCNBS$Cmax
            
        theRow  = p.System_ISCCNBS[ mask.H & mask.V & mask.C,  , drop=FALSE]
        
        if( nrow(theRow) != 1 )
            {
            log_level( ERROR, "(internal).  Expected exactly 1 match for hvc=%g,%g,%g, but found %d.\n", 
                                    hvc[1],hvc[2],hvc[3], nrow(theRow) )
            next
            }
            
        number[i]   = theRow$Number
        name[i]     = munsellinterpol::CentroidsISCCNBS$Name[ theRow$Number ]           # we use the shortcut that rownumber is actually ISCC-NBS Number !
        centroid[i] = munsellinterpol::CentroidsISCCNBS$MunsellSpec[ theRow$Number ]    # we use the shortcut that rownumber is actually ISCC-NBS Number !
        }
        
    out = data.frame( row.names=rnames )
    out$HVC         = HVC
    out$Number      = number
    out$Name        = name
    out$Centroid    = centroid
    
    return( out )
    }
    
    