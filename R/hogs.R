

hogs <- function( iPos=1 )
    {
    theNames    = ls(iPos)
    
    theList = sapply( theNames, function(x) get(x, pos=iPos )  )

    n = length(theList)    
    
    if( n == 0 ) { return(NULL) }    
    
    class1  <- function( x )    { return( class(x)[1] ) }
    
    out     = data.frame( name=theNames, size=sapply( theList, object.size ), 
                                mode=sapply(theList,mode),  class=sapply(theList,class1),  
                                stringsAsFactors=F )

    perm    = order( out$size )
    
    out     = out[ perm, ]
    
    out     = rbind( out, data.frame(name="Total:", size=sum(out$size), mode=NA, class=NA ) )
    
    #   row.names( out )    = theNames
    
    #z[ n+1 ] = sum(z)
    
    #names(z)[n+1] = "Total:"
    
    return( out )
    }
    
removeFunctions  <-  function( iPos=1, .exceptions=c("removeFunctions","hogs") )
    {
    theHogs = hogs()
    
    if( is.null(theHogs) )  return(FALSE)
    
    #   print( theHogs )
     
    df_sub  = subset( theHogs, mode=='function' )
    
    df_sub  = df_sub[ order(df_sub$name), ]
    
    #   print( df_sub )
    
    idx = match( .exceptions, df_sub$name )
    if( 0 < length(idx) )
        #   do not remove these
        df_sub = df_sub[ -idx, ]
        
    #	print( df_sub )

    n   = nrow( df_sub )
    
    if( n == 0 )
        {
        cat( "No functions to remove !\n" )
        return(FALSE)
        }
        
    mess    = sprintf( "About to remove %d functions...\n", n )
    cat( mess )
    
    print( df_sub$name )
    
    keydown <- function(key) { return( key )}

    
    key = readline( prompt="Proceed with removal ?  [Y for Yes]" )
    
    #   print(key)
    
    ok  = toupper(key) == "Y"
    
    if( ! ok )  return(FALSE)
    
    rm( list=df_sub$name, pos=iPos )
    
    return(TRUE)
    }
    