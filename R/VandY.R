

#   YfromV()
#
#   V       a numeric vector of Munsell Values
#           Each V should be in the interval [0,11]
#
#   return  vector of computed Ys, with NA wherever V is invalid
#           this is absolute reflectance Y, as a percentage
#           When V=10, Y=100 exactly.  This is the ASTM-D1535 convention.
#           When V=11, Y=128.7129 approximately
#
#   author: Glenn Davis

YfromV <- function( V, which='ASTM' )
    {
    if( ! missing(which) )
        {
        w  = pmatchYV( which )
        
        if( is.na(w) )
            {
            log.string( WARN, "which='%s' is invalid.", which )
            return( rep(NA_real_,length(V) ) )
            }    
        which   = w
        }
        
    if( which == 'ASTM' )
        out = ( ((((81939*V - 2048400)*V  + 23352000)*V - 22533000)*V + 119140000)*V  ) / 1.e8
    else if( which=='OSA'  ||  which=='MGO' )
        {
        out = ((((8404*V - 210090)*V  + 2395100)*V  - 2311100)*V + 12219000)*V
        
        if( which == 'OSA' )
            out = out / 10256800
        else 
            out = out / 1.e7
        }    
    else if( which == "MUNSELL" )
        {
        b2  = 1474^2
        ac4 = pmin( 4*4740*V^2, b2 )
        out =  50 * ( (1474 - sqrt(b2 - ac4)) / 474 )
        }
    else if( which == "PRIEST" )
        out = V^2
    else
        {
        log.string( FATAL, "Internal Error. which='%s' is invalid.", as.character(which) )
        return(NULL)
        }

    return( out )
    }
    
    
#   VfromY()
#
#   Y   numeric vector of reflectances   
#
#   make a spline function with a large number of 'knots'
# 
VfromY  <- function( Y, which='ASTM' )
    {
    if( ! missing(which) )
        {
        w  = pmatchYV( which )
        if( is.na(w) )
            {
            log.string( WARN, "which='%s' is invalid.", which )
            return( rep(NA_real_,length(Y) ) )
            }    
        which   = w
        }    
        

    if( which %in% names(p.VfromY) )
        out = p.VfromY[[which]](Y)
    else if( which == "MUNSELL" )
        out = sqrt( 1.474*Y - 0.00474*Y^2 )
    else if( which == "PRIEST" )
        out = sqrt(Y)
    else
        {
        log.string( FATAL, "Internal Error. which='%s' is invalid.", as.character(which) )
        return(NULL)
        }        
        
    return( out )
    }
    
    
        
makeVfromYs <- function()
    {
    whichvec = c( 'ASTM', 'OSA', 'MGO' ) 
    
    #log.string( INFO, "Making %d splinefuns...", length(whichvec) )
    #time_start  = gettime()
    
    out = list()
    
    #   these lookup Vs derived by some experimentation - see test-VandY.R for the number of digits of accuracy
    V1  = 0.025 * (-8:119)   # note that 0 is in V1 (this is important so that 0 -> 0).    Old sequence was seq(-0.2,2.98,len=121)
    V2  = seq( 3^(1/2), 10.5^(1/2), len=181 ) ^ (2)
    V   = sort( c( V1, V2, 10 ) )   # ensure that 10 is in V
    
    for( w in whichvec )
        {
        #mess    = sprintf("makeVfromYs().  DEBUG.  Making p.VfromY() for '%s'  (%d Values)...\n", w, length(V) )
        #cat( mess, file=stderr() )
        out[[w]]    = splinefun( YfromV(V,which=w), V, method='fmm' )       # declared in events.R
        }
    
    #log.string( INFO, "done.  [in %g sec]\n", gettime()-time_start )   # less than 0.25 seconds
    
    return( out )
    }
    
    
    
pmatchYV <- function( which )
    {
    full    = c( 'ASTM', 'OSA', 'MGO', 'MUNSELL', 'PRIEST' )

    idx = pmatch( toupper(which), full )
    
    if( is.na(idx) )    return( NA_character_ )
    
    return( full[idx] )
    }
    
    