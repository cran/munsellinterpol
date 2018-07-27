

#   private global variables.  The initial 'p.' means private

#   This group is stored in sysdata.rda
#   p.LookupList        HVC -> xy lookup table(s).  It must be unlocked because the xy of the neutrals (xyC) may be changed during a session.
#   p.InversionCoeffs   3D table of coefficients for approximating A and B as polynomials of a,b.  It stays locked.
#   p.xyC               4x2 matrix of xy coords.  It stays locked.
#   p.sRGB2XYZ          3x3 matrix. It stays locked.
#   p.XYZ2sRGB          3x3 matrix. It stays locked.
#   p.OptimalHull       named list of convex hulls of optimal colors, indexed by 'C', 'D65', etc.  It stays locked.
#   p.xyz1931           the CIE 2-degree color matching functions, at 5nm step. It stays locked.
#   p.ACDs              illuminants A, C, D50, D55, D65, D75, at 5nm step. It stays locked.
#   p.System_ISCCNBS    the elementary blocks for the ISCC-NBS naming system. It stays locked.
#
#   This group is programmatically created during .onAttach()
#   p.VfromY            list of 3 splinefuns; not safe to store them in sysdata.rda.  It must be unlocked.
#   p.microbenchmark    logical value, whether the package microbenchmark is loaded.  It must be unlocked.
#   p.D65toC_CAT        CAT from Illuminant D65 to C.  Used in sRGB conversion.  It must be unlocked.
#   p.CtoD65_CAT        CAT from Illuminant C to D65.  Used in sRGB conversion.  It must be unlocked.

#   These are build options
#   p.vinterpOverride   if TRUE, change from 'cubic' to 'linear' when Value < 2


p.VfromY            = NULL
p.microbenchmark    = FALSE
p.D65toC_CAT        = NULL
p.CtoD65_CAT        = NULL

p.vinterpOverride   = FALSE

    
.onLoad <- function( libname, pkgname )
    {    
    #   unlockBinding() fails here in .onLoad(), so use .onAttach() instead

    }
    
    
.onAttach <- function( libname, pkgname )
    {
    #print( libname )
    #print( pkgname )
    
    if( FALSE )
    {
    info    = library( help='munsellinterpol' )        #eval(pkgname) 
    info    = format( info )
    mask    = grepl( "^(Version|Author|Built)", info )     #Title
    info    = gsub( "[ ]+", ' ', info[mask] )
    mess    = sprintf( "This is %s", pkgname )
    mess    = paste( c( mess, info ), collapse='.  ' )   #; cat(mess)
    packageStartupMessage( mess )
    }
    

    
    unlockBinding( "p.LookupList", asNamespace('munsellinterpol') )         # asNamespace(pkgname) here generates a NOTE !           
    unlockBinding( "p.VfromY", asNamespace('munsellinterpol') )             # asNamespace(pkgname) here generates a NOTE !
    unlockBinding( "p.microbenchmark", asNamespace('munsellinterpol') )     # asNamespace(pkgname) here generates a NOTE ! 
    #   unlockBinding( "p.ListRGB", asNamespace('munsellinterpol') )            # asNamespace(pkgname) here generates a NOTE !     
    #   unlockBinding( "p.InversionCoeffs", asNamespace('munsellinterpol') )    # asNamespace(pkgname) here generates a NOTE !        
    #   packageStartupMessage( "4 variables unlocked." )    
        

    p.microbenchmark    <<- requireNamespace( 'microbenchmark', quietly=TRUE )  #;  cat( "p.microbenchmark=", p.microbenchmark, '\n' )
    
    p.VfromY            <<- makeVfromYs()    #  this loads the list p.VfromY, and takes less than 0.25 seconds
    
    if( requireNamespace( "spacesXYZ", quietly=TRUE ) )
        {
        unlockBinding( "p.D65toC_CAT", asNamespace('munsellinterpol') )
        unlockBinding( "p.CtoD65_CAT", asNamespace('munsellinterpol') )
    
        white.D65   = c( 0.3127, 0.3290, 1 )    # xy are from the official sRGB standard
        white.C     = c( p.xyC['NBS',], 100 )
        
        white.D65   = spacesXYZ::XYZfromxyY( white.D65 )
        white.C     = spacesXYZ::XYZfromxyY( white.C )
    
        p.D65toC_CAT    <<- spacesXYZ::CAT( white.D65, white.C, method='Bradford' )
        p.CtoD65_CAT    <<- spacesXYZ::CAT( white.C, white.D65, method='Bradford' )        
        }
    else
        {
        packageStartupMessage( "ERROR.  Cannot load package spacesXYZ."  )
        }        
    
    
    }

    
    