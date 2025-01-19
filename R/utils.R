

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
#   A   a non-empty numeric Nx3 matrix, or something that can be converted to be one
#
#   returns such a matrix, or NULL in case of error
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

        #   notice hack to make log_level() print name of parent function, and *NOT* prepareNx3()
        log_level( ERROR, "Argument '%s' must be a non-empty numeric Nx%d matrix (with N>0). %s='%s...'",
                                    Aname, M, Aname, mess, .topcall=sys.call(-1L) )
        return(NULL)
        }

    if( ! is.matrix(A) )
        A = matrix( A, ncol=M, byrow=TRUE )

    return( A )
    }



prepareNx3_old  <-  function( A, M=3 )
    {
    ok  = is.numeric(A)  &&  ((length(A) %% M)==0)  &&  0<length(A)
    if( ! ok )
        {
        #print( "prepareNx3" )
        #print( sys.frames() )
        mess    = substr( as.character(A)[1], 1, 10 )
        #arglist = list( ERROR, "A must be a non-empty numeric Nx3 matrix (with N>0). A='%s...'", mess )
        #do.call( log_level, arglist, envir=parent.frame(n=3) )
        #myfun   = log_level
        #environment(myfun) = parent.frame(3)
        log_level( ERROR, "Argument A must be a non-empty numeric Nx%d matrix (with N>0). A='%s...'", M, mess )
        return(NULL)
        }

    if( is.null(dim(A)) )
        A = matrix( A, ncol=M, byrow=TRUE )
    else if( ncol(A) != M )
        A = t(A)

    return( A )
    }
