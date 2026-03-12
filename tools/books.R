
#   this file has functions related to the Munsell books,
#   and intended to be called from savePrivateDatasets()



#   this function reads 3 files and merges them into a single data.frame
#
#   returns a data.frame with 9 columns:
#       HVC             1753x3 matrix, with unique rows, sorted by H,V,C in that order. Neutrals get H=0 and come first.
#       FergusonName    the Ferguson color name; only defined for the first 5 books, and otherwise NA
#       bead            logical, present in the bead book
#       newstudent      <same>
#       plant           <same>
#       rock            <same>
#       soil            <same>
#       glossy          <same>
#       matte           <same>

readAllBooks <- function()
    {
    if( ! requireNamespace('zonohedra') )
        {
        cat( "Cannot load package zonohedra.\n" )
        return(NULL)
        }
    
    books5  = readBooks5()      #; print( str(books5) )

    glossyHVC   = readGlossyBook()

    matteHVC    = readMatteBook()

    ok  = ! is.null(books5)  &&  ! is.null(glossy)  &&  ! is.null(matte)
    if( ! ok )  return(NULL)


    #   add columns to books5 so all 3 can be rbound
    books5$glossy   = FALSE
    books5$matte    = FALSE

    #   add columns to glossy so all 3 can be rbound
    glossy  = data.frame( row.names=1:nrow(glossyHVC) )
    glossy$HVC          = glossyHVC
    glossy$FergusonName = NA_character_
    glossy$bead         = FALSE
    glossy$newstudent   = FALSE
    glossy$plant        = FALSE
    glossy$rock         = FALSE
    glossy$soil         = FALSE
    glossy$glossy       = TRUE      # the only one we truly know !
    glossy$matte        = FALSE


    #   add columns to matte so all 3 can be rbound
    matte  = data.frame( row.names=1:nrow(matteHVC) )
    matte$HVC           = matteHVC
    matte$FergusonName  = NA_character_
    matte$bead          = FALSE
    matte$newstudent    = FALSE
    matte$plant         = FALSE
    matte$rock          = FALSE
    matte$soil          = FALSE
    matte$glossy        = FALSE
    matte$matte         = TRUE      # the only one we truly know !

    df  = rbind( books5, glossy, matte )

    #   sort by H, V, C, in that order
    perm    = order( df$HVC[ ,1], df$HVC[ ,2], df$HVC[ ,3] )

    df     = df[ perm, ]

    #   now find duplicates and merge the duplicate rows. This is the hard part.

    dupvec  = zonohedra::grpDuplicated( df$HVC, MARGIN=1L )

    if( max(dupvec) == 0 )
        {
        #   there are no duplicates, all HVCs are singletons, so we are done
        return( df )
        }

    n   = length(dupvec)        # this is the number of rows, with duplicates

    #  Build output
    count   = attr( dupvec, "ngroups" )[1]       # total rows in output, with no duplicates

    #   idxlist maps from output index to list of input indexes
    idxlist = vector( count, mode='list' )

    #   i indexes the input data.frame df, and the increment is 1 or more
    #   k indexes the output data.frame out, and the increment is always 1
    i   = 1L
    for( k in 1:count )
        {
        if( dupvec[i] == 0 )
            #   just a singleton
            idxlist[[k]]    = i
        else
            idxlist[[k]]    = which( dupvec == dupvec[i] )

        i   = i + length( idxlist[[k]] )
        }

    #   sanity check, i should be n+1L
    if( i != n+1L )
        {
        mess    = sprintf( "ERR. index i=%d, but expected %d.\n", i, n+1L )
        cat( mess )
        return(NULL)
        }


    # print( dupvec[ 1:100 ] )

    #return( idxlist[ (count-10):count ] )

    #   allocate output variable, just to hold results, all will be overwritten !
    out = df[ 1:count, ]
    rownames(out)       = NULL

    out$HVC[ ]  = NA_real_  #  make sure all is overwritten
    rownames(out$HVC)   = NULL

    #   fill in rows of out, one at a time
    for( k in 1:count )
        {
        idx = idxlist[[k]]  # index into df

        if( length(idx) == 1 )
            #   a singleton, so just copy row from df
            out[k, ]    = df[idx, ]
        else
            {
            #   do real work
            dfrow   = data.frame( row.names=idx[1] )    #   make data.frame with 1 row, and the same columns as df

            dfrow$HVC           = df$HVC[ idx[1], , drop=FALSE]     # the HVC are all the same, so just use the first one, idx[1]

            #   find a FergusonName that is *not* NA
            namevec = df$FergusonName[idx]
            j       = which( ! is.na(namevec) )
            if( 0 < length(j) )
                dfrow$FergusonName  = namevec[ j[1] ]       #   actually length(j) should be either 0 or 1
            else
                dfrow$FergusonName  = NA_character_

            if( TRUE )
                {
                #   the last 7 columns are all logical and computed the same way, so use a loop
                for( j in (ncol(dfrow)+1L):ncol(df) )
                    {
                    thename = names(df)[j]
                    dfrow[[ thename ]] = any( df[[ thename ]][idx] )
                    }
                }
            else
                {
                dfrow$bead          = any( df$bead[idx] )
                dfrow$newstudent    = any( df$newstudent[idx] )
                dfrow$plant         = any( df$plant[idx] )
                dfrow$rock          = any( df$rock[idx] )
                dfrow$soil          = any( df$soil[idx] )
                dfrow$glossy        = any( df$glossy[idx] )
                dfrow$matte         = any( df$matte[idx] )
                }

            out[k, ]    = dfrow     # overwrite row k in output
            }
        }
        
    #   check that the rows of out$HVC are unique
    dupvec  = zonohedra::grpDuplicated( out$HVC, MARGIN=1L )

    if( 0 < max(dupvec) )
        {
        mess    = sprintf( "ERR. readAllBooks(). output HVC rows are not unique; there are %d non-singletons.\n",
                                max(dupvec) )
        cat( mess )
        return(NULL)
        }


    neutral = (out$HVC[ ,3] == 0)

    rownames( out )[ neutral ] = MunsellNameFromHVC( out$HVC[neutral, ], format='f', digits=2 )

    rownames( out )[ ! neutral ] = MunsellNameFromHVC( out$HVC[ !neutral, ] )

    #   check that all rownames are not NA
    if( any( is.na(rownames(out)) ) )
        {
        mess    = sprintf( "ERR. readAllBooks(). %d row names are NA.\n", sum( is.na(rownames(out)) ))
        cat( mess )
        return(NULL)
        }

    return( out )
    }


#   read data on the 5 Ferguson books

readBooks5   <- function( path="../inst/extdata/Supplement1_3.6.2024.csv" )
    {
    df  = read.table( path, sep='\t', row.names=NULL, col.names=c("Munsell","FergusonName","books"), fill=TRUE, na.strings='' )
    if( is.null(df) )
        {
        cat( "Cannot read", path, '\n', file=stderr() )
        return(NULL)
        }


    HVC = HVCfromMunsellName( df$Munsell )

    valid   = ! is.na( HVC[ ,2] )

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



#   read 1269 chromatic colors
#   and then prepend 32 achromatic neutral colors

readMatteBook <- function( path="../inst/extdata/mattebook-Joensuu.txt" )
    {
    if( ! requireNamespace('munsellinterpol') )
        {
        cat( "Cannot load munsellinterpol.\n" )
        return(NULL)
        }

    linevec = readLines( path )

    #   get the number of colors from the last line
    colors  = sub( "^[ ]+([0-9]+).*", "\\1", linevec[length(linevec)] )  #;  print( colors )
    colors  = as.integer(colors)

    if( ! is.finite(colors)  ||  colors==0 )
        {
        cat( "Cannot determine chromatic color count, from last line.\n" )
        cat( linevec[length(linevec)], '\n' )
        return(NULL)
        }

    linevec = substring( linevec, 16 )  #   ignore everything after column 16

    #   find the line ending in "RED"
    mask    = grepl( "RED$", linevec )
    if( all( ! mask ) )
        {
        cat( "Cannot find start of data; i.e. the line ending in 'RED'.\n" )
        return(NULL)
        }

    linevec = linevec[ (which(mask)[1]+1):length(linevec) ]

    out = munsellinterpol::HVCfromMunsellName( linevec )

    #   remove any NAs
    mask    = is.finite( out[ ,1] )
    out     = out[mask, ]

    #   check length
    if( nrow(out) != colors )
        {
        mess    = sprintf( "ERR. Read %d chromatic colors, but expected %d.\n", nrow(out), colors )
        cat( mess )
        return(NULL)
        }

    #   cleanup the row names
    rownames(out)   = munsellinterpol::MunsellNameFromHVC( out )

    #   prepend 32 achromatic neutrals
    neutral             = cbind( 0, seq(1.75,9.5,by=0.25), 0 )
    rownames(neutral)   = munsellinterpol::MunsellNameFromHVC(neutral, format='f', digits=2 )

    out                 = rbind( neutral, out )

    return( out )
    }



#   read 1485 chromatic colors
#   and then prepend 37 achromatic neutral colors

readGlossyBook <- function( path="../inst/extdata/glossybook-Centore.txt" )
    {
    if( ! requireNamespace('munsellinterpol') )
        {
        cat( "Cannot load munsellinterpol.\n" )
        return(NULL)
        }

    df = read.csv( path, quote='', comment.char='#' )

    #   check length
    count_expected  = 1485L
    if( nrow(df) != count_expected )
        {
        mess    = sprintf( "ERR. Read %d chromatic colors, but expected %d.\n", nrow(df), count_expected )
        cat( mess )
        return(NULL)
        }

    out = munsellinterpol::HVCfromMunsellName( df$Name )

    #   check for NAs
    mask    = is.na( out[ ,2] )
    if( any(mask) )
        {
        mess    = sprintf( "ERR. Found %d NAs, out of %d names.\n", sum(mask), length(mask) )
        cat( mess )
        return(NULL)
        }

    #   prepend 37 achromatic neutrals
    neutral             = cbind( 0, seq(0.5,9.5,by=0.25), 0 )
    rownames(neutral)   = munsellinterpol::MunsellNameFromHVC( neutral, format='f', digits=2 )

    out                 = rbind( neutral, out )

    return( out )
    }


readGlossyNotations <- function( path="G:/color/historical/Munsell/SpectralReflectancesOf2007MunsellBookOfColorGlossy.txt" )
    {
    linevec = readLines( path )

    #   find the line starting with "Name,"
    mask    = grepl( "^Name,", linevec )
    if( all( ! mask ) )
        {
        cat( "Cannot find start of data; i.e. the line starting with 'Name,'.\n" )
        return(NULL)
        }

    linevec = linevec[ (which(mask)[1]+1):length(linevec) ]

    out = sub( "^([.0-9A-Z/]+).+", "\\1", linevec )

    return( out )
    }



#   HVC     Nx3 matrix with HVC in the rows. These should be unique and taken from a Munsell book.
#   check   check that there are no collisions
#
#   returns a numeric vector

hashHVC <- function( HVC, check=TRUE )
    {
    out = HVC %*% c( 8, 1, 810)

    dim(out) = NULL

    if( check  &&  length(unique(out)) != nrow(HVC) )
        {
        count   = length(unique(out))
        mess    = sprintf( "ERR. hashHVC().  Found %d collisions = %d - %d.\n",
                                nrow(HVC)-count, nrow(HVC), count )
        cat( mess )
        return(NULL)
        }

    names(out)  = rownames(HVC)

    return( out )
    }

nameddiff <- function( x, y )
    {
    mask = x %in% y

    return( x[ !mask ] )
    }

allcomparisons <- function( booklist )
    {
    n   = length(booklist)

    pair    = arrangements::combinations( n, 2 )

    for( k in 1:nrow(pair) )
        {
        j1  = pair[k,1]
        j2  = pair[k,2]

        book1   = booklist[[ j1 ]]
        book2   = booklist[[ j2 ]]

        if( nrow(book1) < nrow(book2) )
            {
            #   book1 is smaller
            cat( names(booklist)[j1], " - ", names(booklist)[j2],  ":  ", length( nameddiff(hashHVC(book1),hashHVC(book2)) ), '\n' )
            }
        else
            {
            #   book2 is smaller
            cat( names(booklist)[j2], " - ", names(booklist)[j1],  ":  ", length( nameddiff(hashHVC(book2),hashHVC(book1)) ), '\n' )
            }
        }

    return( invisible(TRUE) )
    }


