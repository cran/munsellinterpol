
    
#   put fn() between timestamp and the msg    
layout_mine <- structure(
    function(level, msg, namespace="munsellinterpol",
                                    .logcall = sys.call(), .topcall = sys.call(-1), .topenv = parent.frame())
        {
        fn  = deparse1( .topcall[[1L]] )
        
        paste0( attr(level, 'level'), ' [', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '] ', namespace, "::", fn, '(). ', msg )
        },
    generator = quote(layout_mine())
)



#   stop on FATAL
appender_mine <- structure(
    function(lines)
        {
        cat(lines, file = stderr(), sep = '\n' )
        
        #   test for STOP
        if( any( grepl("^(FATAL)",lines ) )  )
            {
            stop( "Stopping in package 'munsellinterpol', because level is FATAL.", call.=FALSE )
            }
        },
    generator = quote(appender_mine())
    )

