    
    
        

#   returns list with
#       omega.from.lambda()
#       lambda.from.omega()
#       responsivity.from.omega[[ ]]    a list with length = channels.  No longer needed, just take deriv=1 of the next one !
#       integral.from.omega[[ ]]        a list with length = channels
#   omega is in [0,1]

makeReparamFunctionList <- function( .lambda, .illum, .xyz, .mode='equalize', .pnorm=2, .Ymax=1 )
    {
    n   = length(.lambda)
    
    if( length(.illum) != n )
        {
        log_level( ERROR, "length(.illum) = %d is incorrect.", length(.illum) )
        return(NULL)
        }
        
    if( ! all( dim(.xyz) == c(n,3) ) )
        {
        log_level( ERROR, "dim(.xyz) = %d,%d is incorrect.", dim(.xyz)[1], dim(.xyz)[2] )
        return(NULL)
        }

    coredata    = matrix( .illum, n, 3 ) * .xyz
            
    #   remove rows that are all 0s
    mask    =   ( 0 < rowSums(abs(coredata)) )
    coredata    = coredata[ mask, ]
    .lambda     = .lambda[ mask ]
    
    step.wl = unique( diff(.lambda) )       
    if( length(step.wl) != 1 )
        {
        log_level( ERROR, "%d wavelengths are not regular.", length(.lambda) )        
        return(NULL)
        }      

    #   check that coredata is full-rank
    singular    = svd( coredata, nu=0, nv=0 )$d
    
    #   log_level( DEBUG, "SVD time = %g sec", as.double(Sys.time()) - time_svd )
    
    thresh  = max( dim(coredata) ) * singular[1] * 2^(-52)
    rank = sum( thresh < singular )
    if( rank <  min(dim(coredata)) )
        {
        log_level( ERROR, "The responsivity matrix   is rank-deficient (rank=%d < %d).", 
                                 rank, min(dim(coredata)) )        
        return(NULL)
        }    
        
    n           = nrow(coredata)
    channels    = ncol(coredata)
    
    lambda.center   = .lambda
    lambda.min      = lambda.center[1]  
    lambda.max      = lambda.center[n]  
        
    #   wavelengths defining the bins
    lambda.break    = seq( lambda.min - step.wl/2, lambda.max + step.wl/2, len=(n+1) )        
        
    out = list()
           
    if( .mode == 'equalize' )
        {
        if( .pnorm == 2 )
            omega   = sqrt( rowSums(coredata*coredata) )
        else if( .pnorm == 1 )
            omega   = rowSums( abs(coredata) )
        else
            omega   = rowSums( abs(coredata)^.pnorm )^(1/.pnorm)
        
        if( any(omega == 0) )
            {
            log_level( ERROR, "Responder has responsivity=0 at 1 or more wavelengths, which is invalid." )
            return(NULL)
            }
            
        omega   = c( 0, cumsum(omega) ) # n+1 of these.
        omega   = omega / omega[n+1]    # n+1 values from 0 to 1.  Not regular.  #;     print( omega )
        
        #   for v 2.6-* added ties=min to suppress warning from regularize.values()
        out$omega.from.lambda   = splinefun( lambda.break, omega, method='monoH.FC', ties=min )      # hyman  monoH.FC  natural
        out$lambda.from.omega   = splinefun( omega, lambda.break, method='monoH.FC', ties=min )      # hyman  monoH.FC  natural
        }
    else if( .mode == 'linear' )
        {
        omega   = (0:n) / n  # n+1 values in regular steps
        
        out$omega.from.lambda   = function( lambda )    { (lambda - lambda.break[1]) / (lambda.break[n+1] - lambda.break[1]) }
        
        out$lambda.from.omega   = function( omega )     {  (1-omega)*lambda.break[1]  +  omega*lambda.break[n+1] }
        }
    else
        {
        log_level( ERROR, ".mode='%s' is invalid.", .mode )
        return(NULL)
        }
        
    out$lambda.break    = lambda.break
    out$omega           = omega
             

    out$integral.from.omega = list()        
    out$omega.from.integral = list()     
    
    #   out$responsivity.from.omega = list()    
    
    s   = .Ymax / (sum( coredata[ , 2] ) )
    
    #   print( s * colSums(coredata) )
    
    for( j in 1:channels )
        {
        integral    = s * c( 0, cumsum( coredata[ , j ] ) )   #   ; print( range(integral) )
        
        #   for v 2.6-* added ties=min to suppress warning from regularize.values()        
        out$integral.from.omega[[j]] <- splinefun( omega, integral, method='monoH.FC', ties=min )   # hyman  monoH.FC  natural
        out$omega.from.integral[[j]] <- splinefun( integral, omega, method='monoH.FC', ties=min )   # hyman  monoH.FC  natural
        
        #   omega.center    = out$omega.from.lambda( lambda.center )
        #   out$responsivity.from.omega[[j]]    = splinefun( omega.center, coredata[ ,j] )  # a tiny bit of extrapolation
            
        #   the next line does not work -- only 1 new function is created  -- I think maybe it is a bug
        #   out$responsivity.from.omega[[j]]    <- function( om ) { return( out$integral.from.omega[[j]]( om, deriv=1 ) ) }
        #   print( str(out$responsivity.from.omega[[j]]) )
        }        
        
    return( out )
    }
    
    
    