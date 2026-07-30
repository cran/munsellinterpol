

#   put fn() between timestamp and the msg
layout_mine <- structure(
    function( level, msg, namespace="munsellinterpol",
                            .logcall = sys.call(), .topcall = sys.call(-1L), .topenv = parent.frame(), .timestamp = Sys.time() )
        {
        if( is.null( .topcall[[1L]] ) )
            fn  = "<top>"
        else
            fn  = deparse1( .topcall[[1L]] )

        fn  = sub( "^.*::", "", fn )    # remove namespace at the beginning, if any

        #paste0( attr(level, 'level'), ' [', format(.timestamp, "%Y-%m-%d %H:%M:%S"), '] ', namespace, "::", fn, '(). ', msg )

        #   use "%-5s" so WARN and INFO are extended to length 5, like all the others, except SUCCESS which is not used
        sprintf( "%-5s [%s] %s::%s(). %s", attr(level,'level'), format(.timestamp, "%Y-%m-%d %H:%M:%S"), namespace, fn, msg )
        },
    generator = quote(layout_mine())
)



appender_mine <- structure(
    function(lines)
        {
        cat( lines, file=stderr(), sep='\n' )

        #   test for STOP - disabled in v 3.5-0
#       if( any( grepl("^(FATAL)",lines ) )  )
#           stop( "Stopping in package 'munsellinterpol', because level is FATAL.", call.=FALSE )
        },
    generator = quote(appender_mine())
    )

