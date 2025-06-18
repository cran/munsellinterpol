
#   HVC0        Nx3 or 1x3 matrix, or character vector of MunsellNotation
#   HVC1        Nx3 or 1x3 matrix, or character vector of MunsellNotation
#   symmetric   logical, whether to symmetrize the difference metric
#
#   returns     numeric N-vector of pairwise differences

NickersonColorDifference <- function( HVC0, HVC1, symmetric=TRUE, coeffs=c(0.4,6,3) )
    {
    if( is.character(HVC0) )
        HVC0 = HVCfromMunsellName(HVC0)
    else  
        {
        HVC0 = prepareNx3( HVC0 )
        if( is.null(HVC0) )  return(NULL)
        }
    
    if( is.character(HVC1) )
        HVC1 = HVCfromMunsellName(HVC1)
    else  
        {
        HVC1 = prepareNx3( HVC1 )
        if( is.null(HVC1) )  return(NULL)
        }
        
    if( nrow(HVC0)==1  &&  1<nrow(HVC1) )
        #   replicate HVC0
        HVC0    = matrix(HVC0,nrow(HVC1),3,byrow=TRUE)
        
    if( 1<nrow(HVC0)  &&  nrow(HVC1)==1 )
        #   replicate HVC1
        HVC1    = matrix(HVC1,nrow(HVC0),3,byrow=TRUE)
        
        
    if( nrow(HVC0) != nrow(HVC1) )
        {
        log_level( ERROR, "nrow(HVC0) = %d != %d = nrow(HVC1).",
                            nrow(HVC0), nrow(HVC1) )
        return(NULL)
        }
        
    ok  = is.numeric(coeffs)  &&  length(coeffs)==3  &&  all( is.finite(coeffs) )  &&  all( 0 < coeffs )
    if( ! ok )
        {
        log_level( ERROR, "argument coeffs is invalid." )
        return(NULL)
        }
        
        
    delta   = abs( HVC0 - HVC1 )
    
    #   scale the Delta Hue column by Chroma
    if( symmetric )
        s   = pmin( HVC0[ ,3], HVC1[ ,3] )
    else
        s   = HVC0[ ,3]
        
    # correct for Delta Hue wraparound
    dHue    = pmin( delta[ ,1], 100 - delta[ ,1] )
        
    delta[ ,1]  = s * dHue
    
    out = delta %*% coeffs
    
    dim(out)    = NULL
        
    return(out)
    }
    