

#   HVC     Nx3 matrix with HVC in the rows
#   books   comma separated strings that designate searched books
#
#   the function depends on the data.frame  MunsellBooks   which is in munsellinterpol.rda

roundHVC    <- function( HVC, books )
    {
    if( is.character(HVC) )
        #   interpret as Munsell notation
        HVC = HVCfromMunsellName(HVC)
    else
        {
        HVC = prepareNx3(HVC)
        if( is.null(HVC) )  return(NULL)
        }

    #   check validity of books
    bookvec = strsplit( books, '[ ,]+' )[[1]]   #; cat( "bookvec =", bookvec, '\n' )
    if( length(bookvec) == 0 )
        {
        log_level( ERROR, "argument books is invalid." )
        return(NULL)
        }

    #   get the full and valid book names
    #   find the "FergusonName" column of MunsellBooks
    j  = match( "FergusonName", colnames(munsellinterpol::MunsellBooks) )
    if( is.na(j) )
        {
        log_level( FATAL, "Internal Error. column 'FergusonName' cannot be found." )
        return(NULL)
        }

    name_book   = colnames(munsellinterpol::MunsellBooks)[ (j+1):ncol(munsellinterpol::MunsellBooks) ]   #; cat( "name_book =", name_book, '\n' )

    idx = pmatch( tolower(bookvec), tolower(name_book), duplicates.ok=TRUE )    #; cat( "idx =", idx, '\n' )
    if( any( is.na(idx) ) )
        {
        i   = which( is.na(idx) )[1]
        log_level( ERROR, "string '%s' in argument 'books' is invalid; it does not match any book.", bookvec[i] )
        return(NULL)
        }

    #   make mask for the valid chips
    mask    = apply( munsellinterpol::MunsellBooks[ , name_book[idx], drop=FALSE ], 1, any )

    sample_search    = munsellinterpol::MunsellBooks[ mask, ]

    log_level( INFO, "searching %d books and %d samples...", length(idx), nrow(sample_search) )

    HVCrnd              = matrix( NA_real_, nrow(HVC), ncol(HVC) )
    colnames(HVCrnd)    = c('H','V','C')

    FergusonName    = rep( NA_character_, nrow(HVC) )

    #   search for closest sample one at a time
    for( k in 1:nrow(HVC) )
        {
        if( any( is.na(HVC[k, ]) ) )    next    # HVC[k, ] is NA
        
        HVCsample   = matrix( HVC[k, ], nrow=nrow(sample_search), ncol=ncol(HVC), byrow=TRUE )
        
        #   symmetric=FALSE results in fewer hue mismatches
        dist        = NickersonColorDifference(  HVCsample, sample_search$HVC, symmetric=FALSE )
        
        #   add very small tie-breaker,  Euclidean distance in V and C/2 multiplied by very small amount, Hue is ignored
        VCsample    = cbind( HVCsample[ ,2], 0.5*HVCsample[ ,3] )
        VCsearch    = cbind( sample_search$HVC[ ,2], 0.5*sample_search$HVC[ ,3] )
        
        dist        = dist +  1.e-6 * sqrt( rowSums( (VCsample - VCsearch)^2 ) )

        if( INFO <= log_threshold( namespace="munsellinterpol" ) )
            {
            #   log the closest 5
            log_level( INFO, "least 5 distances: %s.", paste( sort(dist)[1:5], collapse=", "  ) )
            }

        i   = which.min( dist )

        HVCrnd[k, ]     = sample_search$HVC[i, ]
        FergusonName[k] = sample_search$FergusonName[i]
        }


    rnames  = rownames(HVC)
    
    if( is.null(rnames) )   rnames =  MunsellNameFromHVC( HVC )

    if( any(is.na(rnames))  ||  anyDuplicated(rnames) )
        #   rnames is no good because some are NA or duplicated !  Use trivial names instead.
        #   this should be rare
        rnames = 1:nrow(HVC)   

    out = data.frame( row.names=rnames )

    out$HVC             = HVC
    colnames(out$HVC)   = c('H','V','C')

    out[[ "ISCC-NBS Name" ]]    = ColorBlockFromMunsell( HVC )$Name

    out$MunsellRounded  = MunsellNameFromHVC( HVCrnd, digits=3 )
    out$FergusonName    = FergusonName

    return(out)
    }
