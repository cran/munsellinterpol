
    
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
        log.string( ERROR, "MunsellName is not a character vector with positive length." )
        return(NULL)
        }

    out = matrix( NA_real_, n, 3 )
    colnames(out)   = c( "H", "V", "C" )
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
 
    if( any(mask) )
        {
        out[ mask, 1]    = 0
        out[ mask, 2]    = as.numeric( sub1[mask] )
        out[ mask, 3]    = 0
        }
        
    #   now look for invalid Hue and Value
    #   we must check for Hue because hue_minor may have been set to NA above
    #   we do not really have to check for negative V and C, because the minus sign is not in the pattern
    valid   = is.finite(out[ ,1])  &  (out[ ,2] <= 10)
    
    out[ ! valid, ]  = NA

    return( out )
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
    
    colnames(out) = c('H','V','C')
    
    return(out)

    # return( HVCfromMunsellName( MunsellName )  )
    }



    
#   HueString   a character vector of length n, e.g.  c( '2.5RP', '3.2B' )
#
#   return value    numeric vector of length n    
#                   numeric Hue is ASTM Hue in the interval (0,100]
#
HueNumberFromString  <-  function( HueString )
    {
    # Make all letters in Munsell string upper case
    HueString   = toupper(HueString)  
    
    # Remove whitespace from input Munsell string
    HueString   = gsub( ' |\t', '', HueString )    

    #-----   chromatic colors   ----------#
    pattern = "^([0-9.]+)(R|YR|Y|GY|G|BG|B|PB|P|RP)$" 

    sub12 = sub( pattern, "\\1 \\2", HueString )
    
    mask    = (nchar(HueString) < nchar(sub12) )     # a matched string gets longer          #mask = grepl( pattern, MunsellSpecString )   #   ; print( mask )
    
    n = length(HueString)

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
    
    return( out )
    }
    
#   Hue       a numeric vector, with ASTM Hues, some may be NA    
HueStringFromNumber <- function( Hue, format='g', digits=2 )
    {
    out = rep( NA_character_, length(Hue) )
    
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
        
    out[mask_good] = paste0( formatC(frac,digits=digits,format=format,width=1), namevec[idx+1] )

    #   out[mask_good] = sprintf( "%g%s", frac, namevec[idx+1] )
    
    return(out)
    }
    
    
#   HVC     a numeric vector of length 3
#           or a numeric matrix with 3 columns, with HVC in the rows
#
#   return value:  a character vector    
    
MunsellNameFromHVC <- function( HVC, format='g', digits=2 )
    {    
    HVC = prepareNx3(HVC)
    if( is.null(HVC) )  return(NULL)
        
    out = rep( NA_character_, nrow(HVC) )        
    
    mask_finite = apply( HVC, 1, function(x) { all(is.finite(x)) } )
    
    if( ! any(mask_finite) )  return(out)
        
    #   digits_given    = ( is.finite(digits)  &&  0<=digits )

    if( format == 'f' )  HVC[ ,3] = round( HVC[ ,3], digits=digits )
        
    #   compute chromatic mask
    mask_chromatic  = (0 < HVC[ ,3] )   # negative chroma will convert to neutral
    mask_achromatic = ! mask_chromatic 
    
    mask_chromatic[ ! mask_finite ]   = FALSE 
    mask_achromatic[ ! mask_finite ]  = FALSE   #; print(mask_achromatic)
    
    #   use formatC() 

    out[mask_chromatic]     = paste0( HueStringFromNumber(HVC[mask_chromatic,1], format=format, digits=digits),
                                        ' ', formatC( HVC[mask_chromatic,2], digits=digits, format=format, width=1),
                                        '/', formatC( HVC[mask_chromatic,3], digits=digits, format=format, width=1) )
                                
    if( any(mask_achromatic) )
        out[ mask_achromatic ]  = paste0( "N ", formatC( HVC[mask_achromatic,2], digits=digits, format=format, width=1), '/' )

    return(out)
    }
    
        