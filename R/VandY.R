


#   YfromV_unclamped()
#
#   V       a numeric vector of Munsell Values
#           Values a little outside [0,10] are allowed, because they are useful mathematically,
#               for spline-interpolation purposes in makeVfromYs()
#
#   return  vector of computed Ys, with NA wherever V is invalid
#           this is absolute reflectance Y, as a percentage
#           When V=10, Y=100 exactly.  This is the ASTM-D1535 convention.
#           When V=11, Y=128.7129 approximately
#
#   This one is *NOT* exported.

YfromV_unclamped <- function( V, which='ASTM' )
    {
    if( ! missing(which) )
        {
        w  = pmatchYV( which )

        if( is.na(w) )
            {
            #   in the next call, change .topcall so WARN message names the parent function
            event_level( WARN, "which='%s' is invalid.", as.character(which),
                        class="invalid_argument", extra=c(value=which), .topcall=sys.call(-1L) )
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
            out = out / 1.e7        # MGO
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
        event_level( FATAL, "Internal Error. which='%s' is invalid.", as.character(which),
                            class="internal_error", .topcall=sys.call(-1L) )
        return(NULL)
        }

    return( out )
    }



#   YfromV()
#
#   V       a numeric vector of Munsell Values
#           Each V should be in the interval [0,10].
#
#   return  vector of computed Ys, with NA wherever V is invalid
#           this is absolute reflectance Y, as a percentage
#           When V=10, Y=100 exactly.  This is the ASTM-D1535 convention.
#           When V=11, Y=128.7129 approximately
#
#           negative output values are clamped
#
#   This one *IS* exported.
#   It simply calls YfromV_unclamped() and then clamps.

YfromV <- function( V, which='ASTM' )
    {
    out = YfromV_unclamped( V, which=which )

    if( is.null(out) )  return(out)

    #   check for negative Ys, and clamp if necessary

    bad = is.finite(out)  &  out < 0
    if( any(bad) )
        {
        event_level( WARN, "%d Y(s), of %d, computed to be < 0 (min %.5f); clamped to 0.",
               sum(bad), length(bad), min(out[bad]),
               class = "munsell_clamp",
               extra = list(Y = out[bad], Value = V[bad], indexes = which(bad)) )
        out[bad] = 0
        }

    names(out)  = names(V)

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
            event_level( WARN, "Argument which='%s' is invalid.", as.character(which),
                                class="invalid_argument", extra=c(value=which) )
            return( rep(NA_real_,length(Y) ) )
            }
        which   = w
        }


    if( which %in% names(p.VfromY) )
        out = p.VfromY[[ which ]](Y)    #   p.VfromY is a list of 3 splinefun's
    else if( which == "MUNSELL" )
        out = sqrt( pmax(1.474*Y - 0.00474*Y^2,0) )     # do not allow sqrt() of negative number
    else if( which == "PRIEST" )
        out = sqrt( pmax(Y,0) )                         # do not allow sqrt() of negative number
    else
        {
        event_level( FATAL, "Internal Error. which='%s' is invalid.", as.character(which), class="internal_error" )
        return(NULL)
        }

    #   check for negative Values, and clamp if necessary

    bad = is.finite(out)  &  out < 0
    if( any(bad) )
        {
        event_level( WARN, "%d Munsell Value(s), of %d, computed to be < 0 (min %.5f); clamped to 0.",
               sum(bad), length(bad), min(out[bad]),
               class = "munsell_clamp",
               extra = list(Value = out[bad], Y = Y[bad], indexes = which(bad)) )
        out[bad] = 0
        }

    names(out)  = names(Y)

    return( out )
    }



pmatchYV <- function( which )
    {
    full    = c( 'ASTM', 'OSA', 'MGO', 'MUNSELL', 'PRIEST' )

    idx = pmatch( toupper(which), full )

    if( is.na(idx) )    return( NA_character_ )

    return( full[idx] )
    }




#   makeVfromYs() is called only once, from .onLoad()
#
#   it returns a list of 3 splinefun's

makeVfromYs <- function()
    {
    whichvec = c( 'ASTM', 'OSA', 'MGO' )

    #log_level( INFO, "Making %d splinefuns...", length(whichvec) )
    #time_start  = gettime()

    out = list()

    #   these lookup Vs derived by some experimentation - see test-VandY.R for the number of digits of accuracy

    #   V1 is the 1st sequence, from -0.2 to 3
    #   note that 0:3 are in V1.    The previous sequence was seq(-0.2,2.98,len=121); 0.2 is not a dyadic rational.
    #   0 is important so that 0 -> 0 exactly.
    #   Also a few negative numbers are include so there is no spline-boundary artifact at 0.
    V1  = 0.025 * (-8:120)

    #   V2 is the 2nd sequence, from 3 to 10.5.  There are no integers here, but see V3
    #   It goes past 10 to 10.5 so there is no spline-boundary artifact at 10.
    V2  = seq( 3^(1/2), 10.5^(1/2), len=181 ) ^ (2)
    V2  = V2[-1]    # drop number that is approximately 3, since 3 is already in V1

    #   V3 is simple - the integers 4:10
    V3  = 4:10

    #   combine and sort
    V   = sort( c( V1, V2, V3 ) )

    #   V now has all integers 0:10 - these are the lookup-table V-planes

    for( w in whichvec )
        {
        #mess    = sprintf("makeVfromYs().  DEBUG.  Making p.VfromY() for '%s'  (%d Values)...\n", w, length(V) )
        #cat( mess, file=stderr() )
        out[[w]]    = splinefun( YfromV_unclamped(V,which=w), V, method='fmm' )
        }

    #log_level( INFO, "done.  [in %g sec]\n", gettime()-time_start )   # less than 0.25 seconds

    return( out )
    }

