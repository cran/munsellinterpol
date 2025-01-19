



#   x	a vector containing the x coordinates of the rectangular data grid.  Must be increasing.
#   y	a vector containing the y coordinates of the rectangular data grid.  Must be increasing.
#   z	a matrix containing the z[i,j] data values for the grid points (x[i],y[j]).
#   x0	vector of x coordinates, at which to interpolate
#   y0	vector of y coordinates, at which to interpolate
#
#   Warning:  If x and y are not regular, then the spline is not C^1

bicubicCardinal  <- function( x, y, z, x0, y0 )
    {    
    #   check validity
    lenx    = length(x)
    if( lenx != nrow(z) )
        {
        log_level( ERROR, "length(x)=%d != %d=nrow(z)", length(x), nrow(z) )
        return(NULL)
        }

    leny    = length(y)        
    if( leny != ncol(z) )
        {
        log_level( ERROR, "length(y)=%d != %d=ncol(z)", length(y), ncol(z) )
        return(NULL)
        }
    
    n   = length(y0)
    if( length(x0) != n )
        {
        log_level( ERROR, "length(x0)=%d != %d=length(y0)", length(x0), n )
        return(NULL)
        }
        
    zout    = rep( NA_real_, n )
    
    i   = findInterval( x0, x, rightmost.closed=FALSE ) #; print(i)
    j   = findInterval( y0, y, rightmost.closed=FALSE ) 
    
    zz  = rep( NA_real_, leny )
    
    for( k in 1:n )
        {
        i0  = i[k]
        if( i0==0  || x[lenx] < x0[k] )    next    #   outside the rectangle
        
        j0  = j[k]
        if( j0==0  || y[leny] < y0[k] )    next    #   outside the rectangle
        
        xexact  = (x0[k] == x[i0])
        yexact  = (y0[k] == y[j0])
        
        if( xexact  &&  yexact )
            {
            #   exactly on a grid point.  If one of the neighbors is NA it does not matter
            zout[k] = z[ i0, j0 ]   #; log_level( INFO, "Exactly on grid point %g,%g", x0[k], y0[k] )
            next
            }
        
        if( xexact )
            {
            #   on a gridline, so it reduces to splineCardinal()
            zout[k] = splineCardinal( y, z[i0, ], y0[k] )
            next
            }

        if( yexact )
            {
            #   on a gridline, so it reduces to splineCardinal()
            zout[k] = splineCardinal( x, z[ ,j0], x0[k] )
            next
            }         

        jseq    =  max(j0-1,1) :  min(j0+2,leny)  # ; print(jseq)

        #for( jj in jseq )
        #    zz[jj]  = splineCardinal( x, z[ ,jj], x0[k] )
        
        res = splineCardinal( x, z[ ,jseq], x0[k] )     #; print(zz)
        if( is.null(res) )  next

        z4  = zz
        z4[ jseq ]  = res
        
        zout[k] = splineCardinal( y, z4, y0[k] )
        }
    
    return( list( x=x0, y=y0, z=zout ) )    
    }
    


#   x	a vector containing the x coordinates of the rectangular data grid.  Must be increasing.
#   y	a vector containing the y coordinates of the rectangular data grid.  Must be increasing.
#   z	a matrix containing the z[i,j] data values for the grid points (x[i],y[j]).
#   x0	vector of x coordinates, at which to interpolate
#   y0	vector of y coordinates, at which to interpolate
#
#   the name comes from the fact that this calls stats::approx()

bilinearApprox  <- function( x, y, z, x0, y0 )
    {    
    #   check validity
    lenx    = length(x)
    if( lenx != nrow(z) )
        {
        log_level( ERROR, "length(x)=%d != %d=nrow(z)", length(x), nrow(z) )
        return(NULL)
        }

    leny    = length(y)        
    if( leny != ncol(z) )
        {
        log_level( ERROR, "length(y)=%d != %d=ncol(z)", length(y), ncol(z) )
        return(NULL)
        }
    
    n   = length(y0)
    if( length(x0) != n )
        {
        log_level( ERROR, "length(x0)=%d != %d=length(y0)", length(x0), n )
        return(NULL)
        }
        
    zout    = rep( NA_real_, n )
    
    i   = findInterval( x0, x, rightmost.closed=FALSE ) #; print(i)
    j   = findInterval( y0, y, rightmost.closed=FALSE ) 

    for( k in 1:n )
        {
        i0  = i[k]
        if( i0==0  || x[lenx] < x0[k] )    next    #   outside the rectangle
        
        j0  = j[k]
        if( j0==0  || y[leny] < y0[k] )    next    #   outside the rectangle
        
        xexact  = (x0[k] == x[i0])
        yexact  = (y0[k] == y[j0])
        
        if( xexact  &&  yexact )
            {
            #   exactly on a grid point.  If one of the neighbors is NA it does not matter
            zout[k] = z[ i0, j0 ]   #; log_level( INFO, "Exactly on grid point %g,%g", x0[k], y0[k] )
            next
            }
        
        if( xexact )
            {
            #   on a gridline, so it reduces to approx()
            zout[k] = stats::approx( y, z[i0, ], y0[k] )$y
            next
            }

        if( yexact )
            {
            #   on a gridline, so it reduces to splineCardinal()
            zout[k] = stats::approx( x, z[ ,j0], x0[k] )$y
            next
            }         

        jseq    =  j0 : min(j0+1,leny)  # ; print(jseq)

        zz  = rep( NA_real_, leny )
            
        for( jj in jseq )
            zz[jj]  = stats::approx( x, z[ ,jj], x0[k] )$y
                
        zout[k] = stats::approx( y, zz, y0[k] )$y
        }
    
    return( list( x=x0, y=y0, z=zout ) )    
    }
    

    
    
    
    
    


#   x       N-vector of points, must be increasing.  If the spline is to be C^1, then x must be regular
#   y       N-vector giving values at x.   y can also be an NxP matrix.
#           y can have NA values and it will still work, as long as the subset actually referenced are not NA
#   xout    M-vector of new points at which to evaluate
#
#   return value:
#       if y is an N-vector, then a M-vector 
#       if y is an NxP matrix, then an MxP matrix
#
#   Warning:  If x is not regular, then the spline is not C^1
#
splineCardinal <- function( x, y, xout )
    {
    C   = matrix( c(-0.5,1,-0.5,0,  1.5,-2.5,0,1,  -1.5,2,0.5,0,  0.5,-0.5,0,0), 4, 4 )
    
    #a   = 0.5      this is the only a that has the linear-accuracy property
    #C   = matrix( c(-a,2*a,-a,0,  2-a,a-3,0,1,   a-2,3-2*a,a,0,  a,-a,0,0 ), 4, 4 )
    
    forcemat    = is.null( dim(y) )
    if( forcemat ) 
        #   make y into a matrix with 1 column
        dim(y)  = c(length(y),1)
        
    n   = length(x)
    
    if( nrow(y) != n )
        {
        log_level( ERROR, "length(x) = %d  !=  %d = nrow(y).", length(x), nrow(y) )
        return(NULL)
        }
        
    if( n == 1 )    return( y ) # special case,  maybe shouldn't waste the time checking

    p   = ncol(y)
    
    m   = length(xout)
    out = array( NA_real_, dim=c(m,p) )
    
    idx = findInterval( xout, x, all.inside=TRUE )
    
    for( k in 1:m )
        {
        j   = idx[k]    #; print(j)
        
        u   = (xout[k] - x[j]) / (x[j+1] - x[j])
        u2  = u*u
        U   = c( u*u2, u2, u, 1 )       #; print(U)

        jseq    = pmin( pmax( (j-1):(j+2), 1 ), n )
        Y   = y[ jseq, ,drop=FALSE ]
        
        if( any(is.na(Y)) ) next
        
        #   fixup the ends - by linear extrapolation
        if( j == 1 )
            Y[1, ]    = 2*y[1, ] - y[2, ]     # y[1] - (y[2]-y[1])
        if( j == n-1 )
            Y[4, ]    = 2*y[n, ] - y[n-1, ]   # y[n] + (y[n] - y[n-1])
        
        #print( U %*% C ) 
        #print( Y )
        
        out[k, ]  = (U %*% C) %*% Y
        }
        
    if( forcemat )  
        #   y was originally a plain vector,
        #   so force out to be a plain vector too
        dim(out) = NULL
    
    return( out )
    }

makeSurface3D <- function( interp='bicubic', steps=32  )
    {
    x   = seq( -2, 2, by=1/steps )
    y   = x
    
    Z       = matrix(0,5,5)
    Z[3,3]  = 1
    
    zout    = matrix( NA_real_, length(x), length(y) )
    
    if( interp == 'bicubic' )
        interpfun   = bicubicCardinal
    else
        interpfun   = bilinearApprox
    
    for( i in 1:length(x) )
        {
        zout[i, ]   = interpfun( -2:2, -2:2, Z, rep(x[i],length(y)), y )$z
        }
        
    return( invisible( list( x=x, y=y, z=zout ) ) )
    }
    
    

plotEllipse <- function( asp=3 )
    {
    theta   = seq( 0, 2*pi, len=21 )
    
    x   = asp * cos( theta )
    y   = sin( theta )
    
    plot( x, y, type='n' )
    abline( h=0, v=0 )
    
    thetafine   = seq( 0, 2*pi, len=121 )
    
    xyfine  = splineCardinal( theta, cbind(x,y), thetafine )
    #   yfine   = splineCardinal( theta, y, thetafine )
    
    lines( xyfine[ ,1], xyfine[ ,2], col='red' )
    
    return(TRUE)
    }
    