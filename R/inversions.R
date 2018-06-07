

makeDataForPrediction <- function( Munsell2xy, value, p.LookupList )
    {
    dfsub  = Munsell2xy[ Munsell2xy$V == value , ]    # &  Munsell2xy$real
    
    if( nrow(dfsub) == 0 )  return(NULL)
    
    #   add the illuminant as an aimpoint, which will also affect the model a little bit
    xyC = p.xyC[ 'NBS', ]
    dfsub   = rbind( data.frame(H=0,V=value,C=0,x=xyC[1],y=xyC[2],real=TRUE), dfsub  )
        
    #   add a problem xy point, for debugging
    # dfsub   = rbind( dfsub, data.frame(H=18.75,V=1.488,C=11.69,  x=0.71812,y=0.28188,real=FALSE) )
        
    V.vector    = attr( p.LookupList, "V.vector" )
    H.vector    = attr( p.LookupList, "H.vector" )        
    
    if( value < 1 )
        {
        #   replace xy at non-reals by p.LookupList xy instead.
        #   in some cases they might be the same.        
        iV  = match( value, V.vector )
        
        C.vector    = attr( p.LookupList[[iV]], "C.vector" )
        
        for( j in 1:nrow(dfsub) )
            {
            if( dfsub$real[j] )   next

            iH  = match( dfsub$H[j], H.vector )
            iC  = match( dfsub$C[j], C.vector )
            dfsub$x[j]   = p.LookupList[[iV]]$x[ iH, iC ]
            dfsub$y[j]   = p.LookupList[[iV]]$y[ iH, iC ]                    
            }

        mask    = is.na(dfsub$x) 
        if( FALSE  &&  any(mask) )
            {
            mess    = sprintf( "makeDataForPrediction() WARN.  for value=%g, %d of the rows have x==NA.",
                                        value, sum(mask) )
            cat( mess,'\n', file=stderr() )
            print( dfsub[ mask, ] )
            }
        }    
        
        
        
    #   add A,B
    tmp = ABfromHC( dfsub$H, dfsub$C )
    dfsub$A    = tmp$A
    dfsub$B    = tmp$B
           
    #   add Y    
    dfsub$Y = YfromV(dfsub$V) 
        
    #   add a,b
    XYZ.C   = xyY2XYZ( c( xyC, 100 ) )   

    xyY = cbind( dfsub$x, dfsub$y, dfsub$Y )
    XYZ = xyY2XYZ( xyY )
    Lab = xyz2lab( XYZ, XYZ.C )    
    
    #   dfsub$L    = Lab[ ,1]
    dfsub$a    = Lab[ ,2]
    dfsub$b    = Lab[ ,3]
    
    #   add offset from xyC
    dfsub$xdelta    = dfsub$x - xyC[1]
    dfsub$ydelta    = dfsub$y - xyC[2]
    
    return( dfsub )
    }
    
#   
#   data    as returned from makeDataForPrediction()
#    
#   return value
#       a data.frame with predicted (estimated) columns for AB, HC, xy, and ab.  And these attributes
#       "coeffs"            2x9 matrix of coefficients
addPredictions <- function( data, warn=TRUE )
    {
    #   dfreal  = dfsub[ dfsub$real, ]  

    #   get value from first row
    value   = data$V[1]
    
    out = data    
    
    if( 1 <= value )
        {
        #   use a,b to predict A,B
        modA = lm( A ~ polym(a,b,degree=3,raw=TRUE) + 0, data=data )    
        modB = lm( B ~ polym(a,b,degree=3,raw=TRUE) + 0, data=data )           
        #modA = lm( A ~ a + b + I(a*b) + I(a^2) + I(b^2) + I(a^3) + I(a^2*b) + I(a*b^2) + I(b^3) + 0, data=data )    
        #modB = lm( B ~ a + b + I(a*b) + I(a^2) + I(b^2) + I(a^3) + I(a^2*b) + I(a*b^2) + I(b^3) + 0, data=data )    
        #print( coef(modA) )
        #print( coef(modB) )     
        }
    else
        {
        #   use xdelta,ydelta to predict A,B  
        modA = lm( A ~ polym(xdelta,ydelta,degree=3,raw=TRUE) + 0, data=data )    
        modB = lm( B ~ polym(xdelta,ydelta,degree=3,raw=TRUE) + 0, data=data )   
        }
        
    coeffs  = rbind( stats::coef(modA), stats::coef(modB) )        
        
    colnames(coeffs)    = gsub( ' ', '', names( stats::coef(modA) ) )      # in this case, space are an annoyance
    rownames(coeffs)    = c('A','B')
    #   print( coeffs )
    attr( out, "coeffs" )   = coeffs

    
    
    
    #   add predicted AB
    out$A.pred  = predict( modA, newdata=data )
    out$B.pred  = predict( modB, newdata=data )  
    
    if( warn  &&  any( is.na(out$A.pred) ) )
        log.string( WARN, "%d of A.pred are NA", sum(is.na(out$A.pred)) )

    #   add predicted HC
    tmp = HCfromAB( out$A.pred, out$B.pred )
    out$H.pred  = tmp$H
    out$C.pred  = tmp$C
    
    #   add predicted xy
    HVC = cbind( out$H.pred, value, out$C.pred )
    
    xyY = MunsellToxyY( HVC, warn=FALSE )$xyY
    out$x.pred  = xyY[ ,1]
    out$y.pred  = xyY[ ,2]
    
    if( warn  &&  any( is.na(out$x.pred) ) )
        log.string( WARN, "%d of x.pred are NA", sum(is.na(out$x.pred)) )
    
    #   add predicted ab
    xyC = p.xyC['NBS', ]    
    XYZ.C   = xyY2XYZ( c( xyC, 100 ) )

    xyY = cbind( out$x.pred, out$y.pred, data$Y )
    XYZ = xyY2XYZ( xyY )
    Lab = xyz2lab( XYZ, XYZ.C )    

    out$a.pred = Lab[ ,2]
    out$b.pred = Lab[ ,3]
    
    mask    = is.na(out$a.pred) 
    if( warn  &&  any(mask) )
        {
        log.string( WARN, "%d of a.pred are NA", sum(mask) )
        print( out[mask, ] )
        }
        
    #   print( str(out) )        
        
    return(out)
    }
    
    

    
#   for each of the 15 values, make models for A and B, as polynomials in a and b    
#   return value:
#       list of 15 2x9 matrices.  [[iV]]['A' or 'B'][ index coefficient ]     # 3D array  15x2x9  [iV]['A' or 'B'][ index coefficient ]   
#
makeInversionCoeffs  <-  function( Munsell2xy, p.LookupList, warn=TRUE )
    {
    time_start  = gettime()
    
    V.vector    = sort( unique(Munsell2xy$V) )
    V.vector    = c( 0, V.vector )    
    
    #   triangular  = 4*(4+1)/2     # 10 terms in *full* polynomial, including intercept. like bowling pins.
    
    out = vector( length(V.vector), mode='list' )       #array( NA_real_, dim=c( length(V.vector), 2, triangular-1 ) )
    names(out)  = as.character(V.vector)
    
    #   AB  = c('A','B')
    #   dimnames(out)  = list( as.character(V.vector), AB, NULL )
    
    #   in this loop, start at 2 to skip Value = 0, for which we have no data
    for( iV in 2:length(V.vector) )
        {
        value   = V.vector[iV]
        
        dfv = makeDataForPrediction( Munsell2xy, value, p.LookupList )

        dfv = addPredictions( dfv, warn=warn )
        
        out[[iV]]   = attr( dfv, 'coeffs' )
        
        if( iV == 2 )
            # shrink a little bit because of problems in the very dark purple area
            out[[iV]]   = 0.8 * out[[iV]]  
        }
        
    #   correct the last name
    #   dimnames(out)[[3]] = names( stats::coef(mod.poly) )
    
    #   data for Value=0 is not available
    #   copy from Value=0.2 to Value=0
    out[[1]]    = out[[2]]
    
    time_elapsed    = gettime() - time_start
    mess    = sprintf( "made inversion coeffs in %g seconds.\n", time_elapsed )
    cat(mess)
    
    return(out)
    }