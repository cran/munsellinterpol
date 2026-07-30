
#   level       logging level
#   fmt         format string, passed to sprintf()
#   ...         extra args passed to sprintf()
#   class       a string (or even multiple strings) that can be used to give the type of error that occured
#   extra       named list appended to the custom condition object,
#               with event data that may be useful to the handler
#   .topcall    call, from which the name of the calling function is derived
#
#   *)  constructs single-line message - the record - using current log_layout() callback
#               which is layout_mine() by default
#
#   calls appropriate function, depending on level
#
#   for levels greater than WARN, event_level() and log_level() are equivalent.

event_level <- function( level, fmt, ..., class=character(), extra=NULL, .topcall=sys.call(-1L) )
    {
    if( logger::WARN < level )
        {
        #   just normal logging, no condition
        return( logger::log_level( level, fmt, ... , .topcall=.topcall ) )
        }


    #====   logging work begin  ====#

    message = logger::formatter_sprintf( fmt, ... )

    #   get the current layout function, and call it to compute the layout record
    layout  = logger::log_layout( namespace="munsellinterpol" )

    if( is.symbol(layout) ) layout  = function_from_symbol( layout )

    if( ! is.function(layout) )
        #   cannot log, give up !
        stop( "In package munsellinterpol. Cannot find the logger layout function." )

    record  = layout( level, message, .topcall=.topcall )


    #====   condition work begin    ====#

    extra   = c( list( level=level, package="munsellinterpol", record=record ), extra )

    if( level < logger::WARN )
        #   level must be ERROR or FATAL
        out = event_error( message, class, extra )
    else
        #   level must be WARN
        out = event_warn( message, class, extra )

    return( invisible(out) )
    }




#   message     a short description of the error
#   class       extra class strings, added to the class of the condition object, which can indicate the type of error
#   extra       list of extra data added to condition object, including extra$level and extra$record
#
#   *) constructs custom condition object, with passed class and extra, including extra$record.
#   *) signals the condition
# If condition is not handled, then:
#   *)  get current logger threshold
#   *)  if threshold >= level, writes the record using current log_appender() callback
#               which is appender_mine() by default
#   *) if a restart that applies to "abort" is available, invokes it.  Otherwise stops with the same condition object.

event_error <- function( message, class, extra )
    {
    #   in the next line, call=.topcall leads to expansion of all the args too !  Which I do not want !
    cond    = structure( c( list( message=message, call=NULL ), extra ),
                        class = c( class, "munsellinterpol_error", "error", "condition")
                        )

    signalCondition(cond)       # exiting handlers leave here

    #   condition was unhandled

    #====   logging work resume   ====#

    #   get current logging threshold, and compare with level
    thresh  = logger::log_threshold( namespace="munsellinterpol" )

    appended    = FALSE

    if( extra$level <= thresh )
        {
        #   OK to append/write the record
        #   get the current appender function, and append the layout record
        appender    = logger::log_appender( namespace="munsellinterpol" )

        if( is.symbol(appender) )   appender  = function_from_symbol( appender )

        if( ! is.function(appender) )
            #   cannot log, give up !
            stop( "In package munsellinterpol. Cannot find the logger appender function." )

        appender( extra$record )
        appended    = TRUE
        }

    #====   logging work end   ====#


    #====   condition work resume    ====#

    ra <- computeRestarts("abort")

    if( appended && length(ra) )
        #   The next line silently stops execution.
        #   The silence is OK because the record has already been appended.
        invokeRestart(ra[[1L]])             # usually the same as invokeRestart("abort")
    else
        #   The next line stops execution, and also reports cond$message,
        #   which is the right thing to do since nothing has been reported yet.
        stop( cond )

    #====   condition work end    ====#


    invisible(NULL)
    }




#   message     a short description of the error
#   class       extra class strings, added to the class of the condition object, which can indicate the type of error
#   extra       list of extra data added to condition object, including extra$record.  extra$level == WARN always.
#
#   *) constructs custom condition object, with passed class and metadata, and the record.
#   *) signals the condition
# If condition was not handled with muffleWarning(), then
#   *)  get current logger threshold
#   *)  if threshold >= WARN, writes the record using current log_appender() callback
#               which is appender_mine() by default

event_warn <- function( message, class, extra )
    {
    #   in the next line, call=.topcall leads to expansion of all the args too !  Which I do not want !
    cond    <- structure( c( list( message=message, call=NULL ), extra ),
                        class = c(class, "munsellinterpol_warning", "warning", "condition")
                        )

    withRestarts(
        {
        signalCondition(cond)                   # handlers get first crack

        #====   logging work resumed    ====#
        #   get current logging threshold, and compare with WARN
        thresh  = logger::log_threshold( namespace="munsellinterpol" )

        if( WARN <= thresh )
            {
            #   OK to write/append the record
            #   get the current appender function, and append the layout record
            appender    = logger::log_appender( namespace="munsellinterpol" )

            if( is.symbol(appender) )   appender = function_from_symbol( appender )

            if( ! is.function(appender) )
                #   cannot log, give up
                stop( "In package munsellinterpol. Cannot find the logger appender function." )

            appender( extra$record )
            }
        #====   logging work end    ====#
        },
        #   establish the restart named "muffleWarning".
        #   this is necessary so that suppressWarnings() works with this function, just like it does with warning().
        muffleWarning = function() { return(NULL) }     # so suppressWarnings() works as it does for warning()
        )

    #====   condition work end    ====#

    invisible(NULL)
    }


#   symb    a symbol
#
#   returns a function, or NULL if not found

function_from_symbol    <- function( symb )
    {
    #   in this call to base::get(), inherits=TRUE is necessary
    fun = tryCatch( base::get( symb, mode="function", inherits=TRUE ), error = function(e) { e } )

    if( is.function(fun) ) return( fun )

    #   base::get() failed
    #   try harder using utils::getAnywhere().  This one takes much longer.
    lst = utils::getAnywhere( symb )

    if( 1 <= length(lst$objs)  &&  is.function(lst$objs[[1L]]) )
        {
        #   the first object found is a function; return it
        logger::log_level( logger::DEBUG, "utils::getAnywhere() found symbol '%s' in '%s'.", lst$name, lst$where[1] )
        return( lst$objs[[1L]] )
        }

    return(NULL)
    }
