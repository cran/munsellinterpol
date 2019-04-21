
#   'zonohedron' is an S3 class
#   It presents the zonohedron as an intersection of halfspaces.
#
#   It is stored as a list with
#       W           the original given nx3 matrix, defining n line segments in R^3  (with 1 endpoint at 0).  Must be rank 3.
#       center      middle gray = white/2
#       nonnegative logical.  All entries of W are non-negative
#       face        a data frame with n*(n-1)/2 rows.  Because of central symmetry, only 1/2 of the faces need to be stored.
#                   idx     2 columns with 2 indices to distinct rows of W
#                   normal  to a face of the zonohedron, the cross-product of the 2 rows, then unitized
#                   center  of the face, a parallelogram, *after* centering the zonohedron
#                   beta    the plane constant, *after* centering the zonohedron.  All are positive.

#   zonohedron() constructor for a zonohedron object
#       W       nx3 matrix with rank 3
#       tol     relative tolerance for degenerate faces
#
#   returns: a list as above
#
zonohedron  <- function( W, tol=1.e-6 )
    {
    ok  = is.numeric(W)  &&  is.matrix(W)  &&  3<=nrow(W)  &&  ncol(W)==3
    if( ! ok )
        {
        log.string( ERROR, "argument W is not an nx3 numeric matrix, with n>=3." )
        return(NULL)
        }
        
    if( any( is.na(W) ) )
        {
        log.string( ERROR, "matrix W is invalid because it has %d entries that are NA.", sum(is.na(W)) )
        return(NULL)
        }
    
    p12 = tcrossprod( W[ ,1,drop=F], W[ ,2,drop=F] )
    p13 = tcrossprod( W[ ,1,drop=F], W[ ,3,drop=F] )
    p23 = tcrossprod( W[ ,2,drop=F], W[ ,3,drop=F] )
    
    d12 = p12 - t(p12)
    d13 = t(p13) - p13
    d23 = p23 - t(p23)
    
    n   = nrow(W)
    
    idx = t( utils::combn(n,2) )   # matrix of pairs.  This is kind of slow, but OK for 81 wavelengths
    
    normal  = cbind( d23[idx], d13[idx], d12[idx] )
    
    #   set rows close to 0 to all NAs
    WW  = cbind( W[ idx[,1], ], W[ idx[,2], ] )
    WW2 = .rowSums( WW*WW, nrow(WW), ncol(WW) )
    
    normal2 = .rowSums( normal*normal, nrow(normal), ncol(normal) )
    bad     = normal2 <= tol*tol * WW2
    if( all(bad) )
        {
        log.string( ERROR, "argument W does not have rank 3, with relative tol=%g.", tol )
        return(NULL)
        }


    if( any(bad) )
        {    
        mess    = sprintf( "%d bad normals, out of %d.\n", sum(bad), length(bad) )
        cat(mess)
        #   log.string( INFO, "%d normals flagged as too small, out of %d.", sum(bad), length(bad) )
        normal[ bad, ]  = NA_real_        
        }
    
    #   now unitize
    normal  = normal / sqrt(normal2)    # sqrt(normal2) is replicated to all 3 columns

    out = list()
    
    out$W   = W
    
    out$center  = 0.5 * .colSums( W, nrow(W), ncol(W) )
    
    out$nonnegative = all( 0 <= W )
    
    out$face        = data.frame( row.names=1:nrow(idx) )
    #   out$face$idx    = idx       not really needed for inside()
    out$face$normal = normal
    
    
    if( FALSE )
        {
        #   cluster the unit normal vectors  [reminds me of something similar at Link]
        res = findRowClusters( normal )
        if( is.null(res) )  return(NULL)
        
        out$clusteridx  = res$clusteridx
        out$cluster     = res$cluster
        }    
    
    #   the next line is the biggest bottleneck
    #   the size of the output is on the order n^3
    functional  = tcrossprod( W, normal ) 
    #   pos = 0<S01 ;     S01[ pos ]  = 1 ;     S01[ !pos ] = 0
    
    #   functional is n x n(n-1)/2
    #   each column corresponds to an i<j pair, from array idx[]
    #   the corresponding entries in that column should be very close to 0
    if( TRUE )
        {
        #time_start = gettime()
        
        across  = 1:nrow(idx)
        idx2    = rbind( cbind(idx[,1],across),  cbind(idx[,2],across) )
        
        #   check   = functional[idx2]  ; print(check)      # all these are very small
        
        #   this next line has the effect of replacing an essentially random vertex of the parallelogram
        #   with the center of the parallelogram.
        #   but since all 4 vertices are in the same constant plane, the improvement in the plane constant is unnecessary
        functional[idx2] = 0
        
        #time_elapsed    = gettime() - time_start
        #cat( "time_elapsed =", time_elapsed, '\n' )    # takes  <1% to 5% of total time 
        }
    
    
    #   calculate the n(n-1)/2 face centers of the unit cube,
    #   but then all translated by -0.5 and therefore centered
    #   sign(0) = 0 exactly, so 2 (or more) of the 'spectrum' coords are 0
    face_center = 0.5 * sign(functional)      #;    print( str(S01) )  exact 0 maps to 0

    center = crossprod( face_center, W )        # of the face, a parallelogram, *after* centering the zonohedron
    
    #   out$face$center = center        not really needed for inside()
    
    #   calculate the n(n-1)/2 plane constants beta
    #   these plane constants are for the centered zonohedron
    out$face$beta   = .rowSums( normal * center, nrow(normal), ncol(normal) )
    
    #   since W has rank 3, 0 is in the interior of the _centered_ zonohedron,
    #   and this implies that all n(n-1)/2 beta's are positive, except the NAs.
    #   verify this
    betamin = min( out$face$beta, na.rm=TRUE )
    if( betamin <= 0 )
        {
        log.string( FATAL, "Internal Error.  min(beta)=%g <= 0.", tol, betamin )
        return(NULL)
        }
        
    class(out)  = c( 'zonohedron', class(out) )
        
    return(out)
    }
    
    
##----------        zonohedron methods    -------------##
    
#   x   a zonohedron object
#   g   an Nx3 matrix, etc.
#
#   value   a dataframe with columns
#           within  logical
#           delta   numeric, non-positive means within

inside.zonohedron <- function( x, g )
    {
    g   = prepareNx3( g, 3 )
    if( is.null(g) )    return(NULL)
    
    #   translate g to the centered zonohedron
    gcentered   = g - matrix( x$center, nrow(g), 3, byrow=TRUE ) #; print(gcentered)
    
    hg  = tcrossprod( x$face$normal, gcentered )    #; print( str(hg) )
    
    distance    = abs(hg) - matrix( x$face$beta, nrow(hg), ncol(hg) )
    
    distance    = base::apply( distance, 2, function(z) {suppressWarnings( max(z,na.rm=TRUE) ) } )  #;   print(distance)
    
    distance[ ! is.finite(distance) ]   = NA_real_
    
    if( x$nonnegative )
        {
        #   special override for black
        black   = apply( g, 1, function(v) { isTRUE(all(v==0)) } )  #; print(black)
        if( any(black) )
            distance[black] = 0
            
        #   special override for white.  Fortunately multiplication by 0.5 and 2 preserves all precision.
        white   = apply( g, 1, function(v) { isTRUE(all(v==2*x$center)) } )  #; print(white)
        if( any(white) )
            distance[white] = 0
        }
        
    out = data.frame( within=(distance<=0) )
    rownames(out)   = rownames(g)
    
    #   out$g           = g
    out$delta   = distance

    return(out)
    }
    
    
    
    
#--------       UseMethod() calls           --------------#    
    
inside <- function( x, g )
    {
    UseMethod("inside")
    }    

    
    
#gettime <- function()
#    {
#    return( microbenchmark::get_nanotime() * 1.e-9 )
#    }
    