

#   returns time in seconds, from an arbitrary origin
gettime <- function()
    {
    if( p.microbenchmark )
        return( microbenchmark::get_nanotime() * 1.e-9 )
    else
        return( as.double( base::Sys.time() ) )
    }


#   projectiveMatrix()
#
#   .matrix     invertible matrix, for example a 3x3 matrix with columns the tristimulus coordinates of RGB primaries
#   .unit       non-zero vector.  For example the tristimulus coordinates of white.
#
#   return      square matrix  B, so that
#               B = matrix  %*%  diag(a)  <=>   each column of B is a multiple of the corresponding column in .matrix
#               B %*% 1  =  .unit.      (1 is the vector of all 1s)
#
#   so for colors, B maps RGB to XYZ
#
#    Another way to write these properties:
#        B %*% I = matrix     up to multiples of the columns
#        B %*% 1  =  .unit
#   So I and 1 are the *standard* projective basis,
#   and .matrix and .unit are a different one

projectiveMatrix  <-  function( .matrix, .unit )
    {
    a   = try( solve( .matrix, .unit ), silent=TRUE )

    if( ! is.numeric(a) ) return(NULL)

    ran = range( abs(a) )   #; print(ran)

    if( ran[1] <= 1.e-6 * ran[2] ) return(NULL)

    return( .matrix  %*%  diag(a) )
    }


#   return list with A and B
ABfromHC <- function( H, C )
    {
    theta   = H * pi/50
    list( A = C*cos(theta), B = C*sin(theta) )
    }

HCfromAB <- function( A, B )
    {
    theta   = atan2( B, A )
    list( H=(theta * 50/pi) %% 100, C = sqrt( A^2 + B^2 ) )
    }

hypot<-function(a, b){
# sqrt(a^2 + b^2) without under/overflow.  Author: Jose Gama **/
# http://www.java2s.com/Tutorial/Java/0120__Development/sqrta2b2withoutunderoverflow.htm
r<-0.0
if (abs(a) > abs(b)) {
         r <- b/a
         r <- abs(a)*sqrt(1+r^2)
      } else if (b != 0) {
         r <- a/b
         r <- abs(b)*sqrt(1+r^2)
      }
r
}

###########     argument processing     ##############
#


#   A   a non-empty numeric NxM matrix, or something that can be converted to be one
#
#   returns such a matrix, or NULL in case of error
#
#   This is intended to check user-supplied A, so there is a lot of checking.
#
prepareNx3  <-  function( A, M=3 )
    {
    ok  = is.numeric(A) &&  0<length(A)  &&  (length(dim(A))<=2)

    ok  = ok  &&  ifelse( is.matrix(A), ncol(A)==M, ((length(A) %% M)==0)  )

    if( ! ok )
        {
        mess    = substr( as.character(A)[1], 1, 20 )

        Aname = deparse(substitute(A))

        #   notice .topcall assignment
        #   which makes the logger layout contain the name of the parent function, and *NOT* prepareNx3()
        event_level( ERROR, "Argument '%s' must be a non-empty numeric Nx%d matrix. %s='%s...'",
                                    Aname, M, Aname, mess, class="invalid_argument", .topcall=sys.call(-1L) )
        return(NULL)
        }

    if( ! is.matrix(A) )
        A = matrix( A, ncol=M, byrow=TRUE )

    return( A )
    }


#   HVC     a non-empty numeric Nx3 matrix, with HVC in the rows
#
#   returns the matrix with variables checked and possibly clamped,
#   or NULL in case of error
#
#   This is intended to check user-supplied HVC matrix.
#
#   in all calls to event_level(), note the .topcall assignment
#   which makes the logger layout contain the name of the parent function, and *NOT* prepareHVC()

prepareHVC  <-  function( HVC )
    {
    ok  = is.numeric(HVC)  &&  is.matrix(HVC)  &&  1<=nrow(HVC)  &&   ncol(HVC)==3

    if( ! ok )
        {
        event_level( ERROR, "Argument HVC is not a numeric non-empty Nx3 matrix.",
                                    class="invalid_argument", .topcall=sys.call(-1L) )
        return(NULL)
        }

    #   check Chroma
    bad = HVC[ ,3] < 0
    bad[ is.na(bad) ]   = FALSE
    
    if( any(bad) )
        {
        event_level( WARN, "%d Munsell Chroma(s) (of %d) are < 0 (min Chroma = %.5f); clamped to 0.",
                           sum(bad), length(bad), min(HVC[bad,3]),
                           class = "munsell_clamp",
                           extra = list(Chroma = HVC[bad,3], indexes=which(bad)), .topcall=sys.call(-1L) )
        HVC[bad,3]  = 0
        }

    #   check Value
    bad = HVC[ ,2] < 0
    bad[ is.na(bad) ]   = FALSE
    
    if( any(bad) )
        {
        event_level( WARN, "%d Munsell Value(s) (of %d) are < 0 (min Value = %.5f); clamped to 0.",
                           sum(bad), length(bad), min(HVC[bad,2]),
                           class = "munsell_clamp",
                           extra = list(Value = HVC[bad,2], indexes=which(bad)), .topcall=sys.call(-1L) )
        HVC[bad,2]  = 0
        }

    bad = 10 < HVC[ ,2]
    bad[ is.na(bad) ]   = FALSE
    
    if( any(bad) )
        {
        event_level( WARN, "%d Munsell Value(s) (of %d) are > 10 (max Value = %.5f); clamped to 10.",
                           sum(bad), length(bad), max(HVC[bad,2]),
                           class = "munsell_clamp",
                           extra = list(Value = HVC[bad,2], indexes=which(bad)), .topcall=sys.call(-1L) )
        HVC[bad,2]  = 10
        }

    #   wrap Hue, no clamping necessary
    hue         = HVC[ ,1] %% 100

    # bump 0 to 100 when 0<Chroma, but when Chroma==0, leave Hue==0
    mask        = (hue==0)  &  (0<HVC[ ,3])
    hue[mask]   = 100
    HVC[ , 1]   = hue

    return( HVC )
    }



#   varname name to search for, case-insensitive
#   ...     list of variables
#
#   if varname is found in ..., returns its value
#
#   if varname is not found in ..., returns NULL
#   matching is case-insensitive

intercept   <- function( varname, ... )
    {
    theList =  list( ... )  #; print( theList )

    n   = length(theList)
    if( n == 0 )
        {
        # log_level( ERROR, "No arguments." )
        return(NULL)
        }

    # logger::log_level( logger::DEBUG, "Found %d objects in '...'", n, namespace='munsellinterpol' )

    theNames    = names(theList)

    idx = match( tolower(varname), tolower(theNames) )

    if( is.na(idx) )    return(NULL)

    return( theList[[idx]] )
    }



