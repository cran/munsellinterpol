

#   These datasets are exported and can be inspected by the user.
#   Each one here must be documented.
#   They are all put in the same .RDA file
saveDatasets  <- function( pathout="../data/munsellinterpol.rda" )
    {
    savevec = character(0)

    Munsell2xy  = makeMunsell2xy()
    #   print( str(Munsell2xy) )

    attr( Munsell2xy, "comment" )   = readLines( "../inst/extdata/ReadMe.txt" )

    savevec = c( savevec, "Munsell2xy" )

    pathin    = "../inst/extdata/Centroids_ISCC-NBS.txt"
    CentroidsISCCNBS = read.table( pathin, sep='\t', header=T, stringsAsFactors=F )
    attr( CentroidsISCCNBS, "comment" )  = readLines( pathin, n=23 )

    savevec = c( savevec, "CentroidsISCCNBS" )

    ##  ready to save it
    save( list=savevec, file=pathout, compress='xz' )   #     'xz'  'gzip'  FALSE

    return(TRUE)
    }




#   an advantage of the private data in "sysdata.rda" is that these
#   do not have to be exposed, and therefore documented
savePrivateDatasets  <- function( .path="sysdata.rda" )
    {
    library( logger )
    
    log_threshold( INFO, namespace="munsellinterpol" )
    
    savevec = character(0)

    #   load( "../data/munsellinterpol.rda" )
    #   make little matrix
    p.xyC   = matrix( c( 0.31012,0.31631,  0.3101,0.3163, 0.310,0.316, 0.31006,0.31616 ),  4, 2,    byrow=T )
    # mat = cbind( mat, rep(100,4) )
    rownames(p.xyC) = c( 'JOSA', 'NBS', 'NTSC', 'CIE' )
    colnames(p.xyC) = c('x','y')
    savevec = c( savevec, "p.xyC" )

    if( FALSE )
    {
    #   these two 3x3 matrices are now built in to spacesXYZ
    if( ! requireNamespace('spacesXYZ') )   return(FALSE)

    primary     = matrix( c(0.64,0.33,  0.3,0.6, 0.15,0.06 ), 3, 2, byrow=T )   # from sRGB standard
    primary     = cbind( primary, 1 - rowSums(primary) )
    whiteXYZ    = as.numeric( spacesXYZ::XYZfromxyY( c(0.3127,0.3290,1) ) )       # from sRGB standard
    p.sRGB2XYZ  = projectiveMatrix( t(primary), whiteXYZ )        #; print( p.sRGB2XYZ )
    savevec = c( savevec, "p.sRGB2XYZ" )

    p.XYZ2sRGB  = solve(p.sRGB2XYZ)                             #;  print( p.XYZ2sRGB )
    savevec = c( savevec, "p.XYZ2sRGB" )
    }

    #   Illuminants A, C, D50 D55 D65 D75
    p.ACDs  = read.table( "../inst/extdata/ACDs.5nm.txt", header=TRUE  )
    savevec = c( savevec, "p.ACDs" )

    p.xyz1931   = read.table( "../inst/extdata/xyz1931.5nm.txt", header=TRUE  )
    savevec = c( savevec, "p.xyz1931" )

    p.OptimalHull       = list()

    for( iname in c("C","D65") )
        {
        ispec   = p.ACDs[[ iname ]]
        W   = ispec * as.matrix( p.xyz1931[ ,2:4] )
        W   = 100 * W / colSums(W)[2]
        p.OptimalHull[[ iname ]]  = zonohedron( W )
        }
    savevec = c( savevec, "p.OptimalHull" )



    pathin  = "../inst/extdata/System_ISCC-NBS.txt"
    p.System_ISCCNBS    = read.table( pathin, header=TRUE, sep='\t', stringsAsFactors=F )
    attr( p.System_ISCCNBS, "header" )  = readLines( pathin, n=50 )
    savevec = c( savevec, "p.System_ISCCNBS" )
    
    #   the Munsell books, soil, rock, etc.
    p.Books = readBooks()
    if( is.null(p.Books) )   return(FALSE)
    savevec = c( savevec, "p.Books" )

    if( TRUE )
    {
    Munsell2xy  = makeMunsell2xy()          # this one does not generate a warning

    p.LookupList = munsellinterpol:::makeLookupList( Munsell2xy, p.xyC['NBS', ], kfactor=c(0.70,0.50) )     # default is kfactor = c(0.9,0.5)
    savevec = c( savevec, "p.LookupList" )    # this compresses *very* well

    p.InversionCoeffs =  munsellinterpol:::makeInversionCoeffs( Munsell2xy, p.LookupList, warn=FALSE )
    savevec = c( savevec, "p.InversionCoeffs" )
    }

    ##  finally ready to save it
    save( list=savevec, file=.path, compress='xz' )   #     'xz'  'gzip'  FALSE

    mess    = sprintf( "Saved the following to '%s'.", .path )
    cat( mess, '\n', file=stderr() )
    cat( savevec, '\n', file=stderr() )

    return( invisible(TRUE) )
    }

pingDatasets  <- function( .path="sysdata.rda", .verbose=FALSE )
    {
    theName     = load(.path,verbose=.verbose)
    print( theName )

    if( 0 < length(theName)  &&  .verbose )
        {
        for( k in 1:length(theName ) )
            {
            obj = get( theName[k] )
            cat( '\n', theName[k], '\n' )
            print( str(obj) )
            }
        }

    return( invisible(T) )
    }


makeMunsell2xy <- function(  )
    {
    time_start  = munsellinterpol:::gettime()

    realMunsell = read.table( "../inst/extdata/real.dat", header=T, sep='', stringsAsFactors=F )    # ; print( str(realMunsell) )
    colnames(realMunsell)[1] = 'H'      # was 'h'
    realMunsell$H = HueNumberFromString( realMunsell$H )

    verydark = read.table( "../inst/extdata/verydark.dat", header=T, sep='', stringsAsFactors=F )   # ; print( str(verydark) )
    verydark$H = HueNumberFromString( verydark$H )

    allMunsell  = read.table( "../inst/extdata/all.dat", header=T, sep='', stringsAsFactors=F )     # ; print( str(allMunsell) )
    allMunsell$H = HueNumberFromString( allMunsell$H )

    #   diffMunsell( realMunsell, allMunsell )

    real    = contained( allMunsell, realMunsell )  |  contained( allMunsell, verydark )    # slow, not optimized

    allMunsell$C = as.numeric( allMunsell$C )   # do not want integer here

    out = allMunsell[ -6 ]  # drop Y

    out  = cbind( out, real=real  )

    #   make slight tweak(s)
    idx = which( out$H==72.5  &  out$V==10  &  out$C==2 )
    if( length(idx) == 1 )
        {
        out$x[idx]  = 0.2970
        out$y[idx]  = 0.3075
        }



    mess    = sprintf( "makeMunsell2xy(). INFO.  Made Munsell2xy in %g seconds.", munsellinterpol:::gettime() - time_start )
    cat( mess, '\n', file=stderr() )

    return( invisible(out) )
    }


contained <- function( .data1, .data2 )
    {
    #   brute force compare, not optimized
    n   = nrow(.data1)

    out     = logical(n)

    name1   = deparse( substitute( .data1 )  )
    name2   = deparse( substitute( .data2 )  )

    for( k in 1:n )
        {
        row1    = .data1[ k, ]

        mask    = .data2$H==row1$H  &  .data2$V==row1$V  &  .data2$C==row1$C
        row2    = .data2[ mask, ]

        if( nrow(row2) == 0 )   next

        if( 2 <= nrow(row2) )
            {
            mess    = sprintf( "'%s' is invalid. Duplicated data\n", name2 )
            cat( mess )
            print( row2 )
            return(NULL)
            }

        #   exactly 1 match
        out[k]  = TRUE
        }

    return(out)
    }
    
    
    


readBooks   <- function( path="../inst/extdata/Supplement1_3.6.2024.csv" )
    {
    df  = read.table( path, sep='\t', row.names=NULL, col.names=c("Munsell","color.name","books"), fill=TRUE, na.strings='' )
    if( is.null(df) )
        {
        cat( "Cannot read", path, '\n', file=stderr() )
        return(NULL)
        }
        
    #print( df[1:20, ] )
        
    HVC = HVCfromMunsellName( df$Munsell )
    
    valid   = ! is.na( HVC[ ,1] )
    
    #   take out the non-valid rows
    row.names   = df$Munsell[valid]
    
    df  = df[ valid, ] 
    rownames(df)    = row.names
    
    #   replace the first column with the matrix HVC
    df[[1]] = HVC[ valid, ]
    colnames(df)[1] = "HVC"
    
    #print( df[1:20, ] ) ; print( str(df) )
    
    #   replace the column 'books' with 5 distinct logical columns
    book    = c( B="bead", N="newstudent", P="plant", R="rock", S="soil" )
    
    for( nam in names(book) )
        {
        df[[ book[nam] ]] = grepl( nam, df$books )
        }
        
    #   remove the 'books' column, no longer needed
    df$books    = NULL
        
    #print( df[1:20, ] ) ; print( str(df) )
    
    return( invisible(df) )
    }
        
    
    
    
    
######################      deadwood below      #################################    



diffMunsell <- function( .data1, .data2 )
    {
    #   brute force compare, not optimized
    out     = NULL
    count   = 0

    name1   = deparse( substitute( .data1 )  )
    name2   = deparse( substitute( .data2 )  )

    mat1    = matrix( 0, 1, 2 )
    colnames(mat1)  = c( 'x', 'y' )
    class(mat1) = "model.matrix"
    mat2    = mat1

    for( k in 1:nrow(.data1) )
        {
        row1    = .data1[ k, ]

        mask    = .data2$H==row1$H  &  .data2$V==row1$V  &  .data2$C==row1$C
        row2    = .data2[ mask, ]

        if( nrow(row2) == 0 )   next

        if( 2 <= nrow(row2) )
            {
            mess    = sprintf( "'%s' is invalid. Duplicated data\n", name2 )
            cat( mess )
            print( row2 )
            return(NULL)
            }

        #   exactly 1 match
        count   = count + 1

        mat1[ 1, 1 ]    = row1$x
        mat1[ 1, 2 ]    = row1$y
        mat2[ 1, 1 ]    = row2$x
        mat2[ 1, 2 ]    = row2$y

        ok  = all( mat1 == mat2 )

        if( ! ok )
            {
            delta   = mat1 - mat2
            #   colnames(delta) = c( "x", "y" )      redundant
            class(delta) = "model.matrix"
            out = rbind( out, data.frame( H=row1$H, V=row1$V, C=row1$C, xy1=mat1, xy2=mat2, delta=delta, stringsAsFactors=F ) )
            }
        }


    colnames( out )[4] = name1
    colnames( out )[5] = name2


    mess    = sprintf( "Chips that %s and %s have in common = %d\n", name1, name2, count )
    cat(mess)
    mess    = sprintf( "Chips that are in %s, but not in %s = %d\n", name1, name2, nrow(.data1) - count )
    cat(mess)
    mess    = sprintf( "Chips that are in %s, but not in %s = %d\n", name2, name1, nrow(.data2) - count )
    cat(mess)

    return( out )
    }
