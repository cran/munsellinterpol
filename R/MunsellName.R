

#   HVCfromMunsellName()
#
#   convert Munsell notation to numeric Hue,Value,Chroma
#
#   arguments:
#       MunsellName   a character vector of length N
#
#   return value: an Nx3 matrix with the computed HVCs placed in the rows
#       Hue     is ASTM Hue in the interval (0,100]  (except Hue=0 for neutrals)
#       Value   is in the interval [0,10]
#       Chroma  is positive  (except Chroma=0 for neutrals)
#
#   rownames(out) are set to MunsellName
#
#   author:  Glenn Davis

HVCfromMunsellName <- function( MunsellName )
    {
    n = length(MunsellName)

    if( n==0  ||  ! is.character(MunsellName) )
        {
        event_level( ERROR, "Argument 'MunsellName' is invalid. It must be character vector with positive length.",
                                    class="invalid_argument", extra=list(MunsellName=MunsellName) )
        return(NULL)
        }

    out = matrix( NA_real_, n, 3 )
    colnames(out)   = c( "H", "V", "C" )

    #thenames    = names(MunsellName)
    #if( is.null(thenames) ) thenames    = MunsellName

    rownames(out)   = MunsellName

    # Remove whitespace from input Munsell string
    MunsellName <- gsub( ' |\t', '', MunsellName )

    # Make all letters in Munsell string upper case
    MunsellName <- toupper(MunsellName)

    #-----   chromatic colors   ----------#
    pattern = "^([0-9.]+)(R|YR|Y|GY|G|BG|B|PB|P|RP)([0-9.]+)/([0-9.]+)$"

    #   time_start  = as.double( Sys.time() )

    sub1234 = sub( pattern, "\\1 \\2 \\3 \\4", MunsellName )

    mask    = (nchar(MunsellName) < nchar(sub1234) )     # a matched string gets longer          #mask = grepl( pattern, MunsellName )   #   ; print( mask )

    #   change NAs to FALSE
    mask[ is.na(mask) ] = FALSE

    if( any(mask) )
        {
        #   make a little LUT
        hue = seq( 0, 90, by=10 )
        names(hue)  = c("R","YR","Y","GY","G","BG","B","PB","P","RP")

        sub1234    = sub1234[mask]

        dat     = unlist( strsplit( sub1234, ' ', fixed=T ) ) #;  print( dat )
        dat     = matrix( dat, length(sub1234), 4, byrow=T ) #;    print( dat )

        hue_minor   =  as.numeric(dat[ ,1])
        valid       =  (hue_minor <= 10)
        hue_minor[ ! valid ]    = NA

        out[ mask, 1]   = hue[ dat[ ,2] ]  +  hue_minor
        out[ mask, 2]   = as.numeric( dat[ ,3] )
        out[ mask, 3]   = as.numeric( dat[ ,4] )
        }

    #   cat( "Elapsed: ", as.double( Sys.time() ) - time_start, '\n' )

    #------   achromatic colors     ----#
    pattern = "^N([0-9.]+)(/|/0)?$"

    sub1    = sub( pattern, "\\1", MunsellName )

    mask    = (nchar(sub1) < nchar(MunsellName) ) # a matched string gets shorter.    #mask = grepl( pattern, MunsellName )   #   ; print( mask )

    #   change NAs to FALSE
    mask[ is.na(mask) ] = FALSE

    if( any(mask) )
        {
        out[ mask, 1]    = 0
        out[ mask, 2]    = as.numeric( sub1[mask] )
        out[ mask, 3]    = 0
        }

    #   now look for invalid Hue and Value
    #   we must check for Hue because hue_minor may have been set to NA above
    #   we do not really have to check for negative V and C, because the minus sign is not in the pattern
    bad =  !( is.finite(out[ ,1])  &  (out[ ,2] <= 10) )

    if( any(bad) )
        {
        #   set all relevant rows to NA
        out[ bad, ] = NA

        MunsellName_org = rownames(out) # given original names

        event_level( WARN, "%d Munsell names, out of %d, could not be converted to HVC.", sum(bad), length(bad),
                        class="incomplete", extra=list( invalid=MunsellName_org[bad], indexes=which(bad) ) )
        }

    return( out )
    }





#   HueString   a character vector of length n, e.g.  c( '2.5RP', '3.2B' )
#
#   return value    numeric vector of length n
#                   numeric Hue is ASTM Hue in the interval (0,100]
#
HueNumberFromString  <-  function( HueString )
    {
    n = length(HueString)

    if( n==0  ||  ! is.character(HueString) )
        {
        event_level( ERROR, "Argument 'HueString' is invalid. It must be character vector with positive length.",
                                    class="invalid_argument", extra=list(HueString=HueString) )
        return(NULL)
        }

    HueString_org   = HueString     # save originals for potential WARN later

    # Make all letters in Munsell string upper case
    HueString   = toupper(HueString)

    # Remove whitespace from input Munsell string
    HueString   = gsub( ' |\t', '', HueString )

    #-----   chromatic colors   ----------#
    pattern = "^([0-9.]+)(R|YR|Y|GY|G|BG|B|PB|P|RP)$"

    sub12 = sub( pattern, "\\1 \\2", HueString )

    mask    = (nchar(HueString) < nchar(sub12) )     # a matched string gets longer          #mask = grepl( pattern, MunsellSpecString )   #   ; print( mask )

    out = rep( NA_real_, n )

    if( any(mask) )
        {
        #   make a little LUT
        hue = seq( 0, 90, by=10 )
        names(hue)  = c("R","YR","Y","GY","G","BG","B","PB","P","RP")

        sub12    = sub12[mask]

        dat     = unlist( strsplit( sub12, ' ', fixed=T ) ) #;  print( dat )
        mat     = matrix( dat, length(sub12), 2, byrow=T ) #;    print( dat )

        hue_minor   =  as.numeric(mat[ ,1])
        valid       =  (hue_minor <= 10)
        hue_minor[ ! valid ]    = NA

        out[ mask ]   = hue[ mat[ ,2] ]  +  hue_minor
        }

    bad =  is.na( out )
    if( any(bad) )
        {
        event_level( WARN, "%d Hue strings, out of %d, could not be converted to numerical Hue.", sum(bad), length(bad),
                        class="incomplete", extra=list( invalid=HueString_org[bad], indexes=which(bad) ) )
        }

    return( out )
    }

#   Hue     a numeric vector, with ASTM Hues, some may be NA_real_
#           Hue values are automatically wrapped to (0,100]
#
#   returns character vector of the same length; NA_real_  converted to NA_character_

HueStringFromNumber <- function( Hue, format='g', digits=2 )
    {
    out = rep( NA_character_, length(Hue) )

    names(out)  = names(Hue)

    mask_good   = is.finite(Hue)

    if( ! any(mask_good) )  return(out)

    #   digits_given    = ( is.finite(digits)  &&  0<=digits )
    #   if( digits_given )  Hue = round( Hue, digits=digits )   not needed

    Hue   = Hue[mask_good] %% 100

    Hue[ Hue == 0 ] = 100

    idx     = as.integer( Hue / 10 )
    frac    = Hue - 10 * idx

    mask    = (frac == 0)
    if( any(mask) )
        {
        idx[mask]   = idx[mask]-1
        frac[mask]  = 10
        }

    namevec  = c("R","YR","Y","GY","G","BG","B","PB","P","RP")

    #   use formatC()
    out[mask_good] = paste0( formatC(frac,digits=digits,format=format,width=1,decimal.mark='.'), namevec[idx+1L] )

    #   out[mask_good] = sprintf( "%g%s", frac, namevec[idx+1] )

    return(out)
    }


#   HVC     Nx3 numeric matrix, with HVC in the rows
#           or a numeric vector of length a multiple of 3
#
#   return value:  a character N-vector

MunsellNameFromHVC <- function( HVC, format='g', digits=2, ctol=0 )
    {
    HVC = prepareNx3(HVC)
    if( is.null(HVC) )  return(NULL)

    #   check format
    ok  = is.character(format)  &&  length(format) %in% 1:3
    if( ! ok )
        {
        event_level( ERROR, "Argument 'format' is invalid.", class="invalid_argument", extra=list(format=format) )
        return( NULL )
        }

    if( length(format) == 1 )
        format      = rep( format, 3 )
    else if( length(format) == 2 )
        format[3]   = format[2]

    #   check digits
    ok  = is.numeric(digits)  &&  length(digits) %in% 1:3
    if( ! ok )
        {
        event_level( ERROR, "Argument 'digits' is invalid.", class="invalid_argument", extra=list(digits=digits) )
        return( NULL )
        }

    if( length(digits) == 1 )
        digits      = rep( digits, 3 )
    else if( length(digits) == 2 )
        digits[3]   = digits[2]


    out = rep( NA_character_, nrow(HVC) )

    #   names(out)  = rownames(HVC)

    mask_finite = apply( HVC, 1, function(x) { all(is.finite(x)) } )

    if( ! any(mask_finite) )
        {
        event_level( WARN, "None of the %d HVC rows could be converted to Munsell names; returning all NAs.", nrow(HVC),
                            class="incomplete", extra=c(invalid=HVC) )
        return(out)
        }

    #   snap Chroma to 0 if appropriate
    Chroma  = HVC[ ,3]

    if( format[3] == 'f' )  Chroma = round( Chroma, digits=digits[3] )

    #   if ctol=0, which is the default, then the next line does nothing,
    #   unless a Chroma is negative, when it is set to 0, which is appropriate
    Chroma[ Chroma < ctol ] = 0

    #   compute chromatic and neutral masks
    mask_chromatic  = 0 < Chroma
    mask_neutral    = ! mask_chromatic

    #   force skipping any row that is not finite
    mask_chromatic[ ! mask_finite ] = FALSE
    mask_neutral[ ! mask_finite ]   = FALSE


    if( any(mask_chromatic) )
        #   use formatC() for Value and Chroma
        out[mask_chromatic] = paste0( HueStringFromNumber(HVC[mask_chromatic,1], format=format[1], digits=digits[1] ),
                                ' ', formatC( HVC[mask_chromatic,2], digits=digits[2], format=format[2], width=1, decimal.mark='.' ),
                                '/', formatC( Chroma[mask_chromatic], digits=digits[3], format=format[3], width=1, decimal.mark='.' ) )

    if( any(mask_neutral) )
        #   only Value, the middle column, is used here
        out[ mask_neutral ]  = paste0( "N ",
                                    formatC( HVC[mask_neutral,2], digits=digits[2], format=format[2], width=1, decimal.mark='.'), '/' )

    bad = is.na( out )
    if( any(bad) )
        event_level( WARN, "%d, of %d, HVC rows could not be converted to Munsell names.", sum(bad), length(bad),
                            class="incomplete", extra = c( invalid=HVC[bad, ,drop=FALSE], indexes=which(bad) ) )

    return(out)
    }




#   for compatibility with package colorscience
#
#   MunsellName     character N-vector of Munsell notations
#
#   return value
#       Nx3 matrix of character strings
#
MunsellHVC  <- function( MunsellName )
    {
    out = strsplit( MunsellName, "[ ]+|/" )

    fixup <- function( vec )
        {
        if( vec[1]=='N'  &&  vec[length(vec)] != '0' )  vec = c( vec, '0' )

        if( length(vec) == 3 )
            return(vec)
        else
            return( rep(NA_character_,3) )
        }

    out = lapply( out, fixup )      #; print(out)

    n   = length(out)

    out = unlist( out )

    #   if( length(out) != 3*n )  out = rep( NA_character_, 3*n )

    out = matrix( out, nrow=n, 3, byrow=TRUE )

    rownames(out)   = names( MunsellName )
    colnames(out)   = c('H','V','C')

    return(out)
    }

