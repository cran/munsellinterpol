
FATAL <- 1L
names(FATAL) <- "FATAL"
ERROR <- 2L
names(ERROR) <- "ERROR"
WARN <- 4L
names(WARN) <- "WARN "
INFO <- 6L
names(INFO) <- "INFO "
DEBUG <- 8L
names(DEBUG) <- "DEBUG"
TRACE <- 9L
names(TRACE) <- "TRACE"

log.string <- function( level, msg, ... )
    {    
    conn    = stderr()
    
    if( ! is.integer(level) )
        {
        cat( "log.string(). level is not an integer.\n" )
        return( invisible(FALSE) )
        }
        
    msg = sprintf( msg[1], ... )    # should this really be msg[1] ?

    #   find the name of the calling function
    #print( "log.string" )
    # print( sys.status() )
        
    where   = sys.parent(1)  # ; print(where)
  
    if( 0 < where )
        namefun = tryCatch( deparse(sys.call(where)[[1L]]), error=function(e) "[console]" )
    else
        namefun = "[??]"

    if( ! grepl( "^munsellinterpol", namefun ) )
        namefun = paste0( "munsellinterpol::", namefun, collapse='' )        
        
    mess    = paste0( namefun, "(). ", names(level), ".  ", msg, collapse='' )

    if( level <= FATAL )
        stop( mess, '\n',  call.=FALSE )

    cat( mess, '\n', file=conn ) ;   flush(conn)

    return( invisible(TRUE) )
    }