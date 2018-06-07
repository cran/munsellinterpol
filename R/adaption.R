

#   xyY.src     xyY in the source viewing environment, a 3-vector or an nx3 matrix
#   white.src   xyY of source white
#   white.dest  xyY of destination white
#
#   returns  XYZ adapted to white.dest
adaptxyY  <-  function( xyY.src, white.src, white.dest, method="bradford" )
    {
    #   convert all 3 to XYZ
    xyY.src = prepareNx3( xyY.src )    
    if( is.null(xyY.src) )  return(NULL)
    
    XYZ.src     = xyY2XYZ( xyY.src )
    white.src   = xyY2XYZ( white.src )
    white.dest  = xyY2XYZ( white.dest )
    
    if( is.null(XYZ.src)  ||  is.null(white.src)  ||  is.null(white.dest) ) return(NULL)
    
    XYZ.dest    = adaptXYZ( XYZ.src, white.src, white.dest, method=method )
    
    if( is.null(XYZ.dest) )  return(NULL)
    
    #   now back to xyY
    xyY.dest = XYZ2xyY( XYZ.dest )
    
    #   do an additional test to correct possible roundoff error
    mask    = (xyY.src[ ,1] == white.src[1])  &  (xyY.src[ ,2] == white.src[2])
    if( any(mask,na.rm=TRUE) )
        {
        #   for some rows, xyY.src has the exact chromaticity of white.src
        #   therefore, xyY.dest must have exact chromaticity of white.dest
        xyY.dest[mask,1] = white.dest[1]
        xyY.dest[mask,2] = white.dest[2]
        }
    
    return( xyY.dest )
    }



#   XYZ.src     XYZ in the source viewing environment, a 3-vector or an nx3 matrix
#   white.src   XYZ of source white
#   white.dest  XYZ of destination white
#
#   returns  XYZ adapted to white.dest
adaptXYZ  <-  function( XYZ.src, white.src, white.dest, method="bradford" )
    {
    XYZ.src = prepareNx3( XYZ.src )    
    if( is.null(XYZ.src) )  return(NULL)

    chad = ChromaticAdaptionMatrix( white.src, white.dest, method )      #;   print( chad )
    
    if( is.null(chad) ) return(NULL)
    
    out = tcrossprod( XYZ.src, chad )
    
    #   for( k in 1:nrow(out) ) { out[ k, ]   = chad %*% as.double( out[k, ] ) }

    rownames(out)   = rownames(XYZ.src)
    colnames(out)   = c('X','Y','Z')
    
    return( out )
    }





#   returns a matrix that maps .white_src to .white_dest    
ChromaticAdaptionMatrix <- function( .white_src, .white_dest, .method="bradford" )
    {
    full    = c( "bradford", "vonkries", "mcat02", "scaling" )
    
    idx     = pmatch( tolower(.method), full )
    if( is.na(idx) )
        {
        log.string( ERROR, "method='%s' is invalid\n", .method )
        return(NULL)
        }    
    .method = full[idx]
    
    if( .method == "bradford" )
        Ma = matrix( c(0.8951,0.2664,-0.1614,  -0.7502,1.7135,0.0367,  0.0389,-0.0685,1.0296), 3, 3, byrow=T )
    else if( .method == "vonkries" )
        Ma = matrix( c(0.40024,0.7076,-0.08081,  -0.2263,1.16532,0.0457,  0,0,0.91822), 3, 3, byrow=T )
    else if( .method == "mcat02" )
        Ma = matrix( c( 0.7328, 0.4296, -0.1624,  -0.7036, 1.6975, 0.0061, 0.0030, 0.0136, 0.9834 ), 3, 3, byrow=T )
    else if( .method == "scaling" )
        Ma = diag(3)

    crm_src  = coneResponseMatrix( Ma, .white_src )
    
    crm_dest = coneResponseMatrix( Ma, .white_dest )
    
    if( is.null(crm_src)  ||  is.null(crm_dest) )   return(NULL)
    
    #   print( lms_dest )
    
    out = solve(crm_dest) %*% crm_src

    return( out )
    }
    
    
#   .Ma     the matrix for the method
#   .white  the reference white    
#   returns a matrix that maps .white to (1,1,1)
coneResponseMatrix  <-  function( .Ma, .white )
    {   
    lms = .Ma %*% as.numeric(.white)
    
    return( diag( as.double(1/lms) ) %*% .Ma )
    }
        