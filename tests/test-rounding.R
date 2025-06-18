
library( munsellinterpol )

options( width=144 )

#
#   checkColorLookup()
#   performs lookup with ColorBlockFromMunsell() and verifies that the color block number is correct
#
#   return value:  TRUE or FALSE
#
testRounding <- function( count=100 )
    {
    set.seed(0)
    
    #   make hues that are in the book
    countbad    = 0
    
    #   make vector of Hues that actually appear in the soil book
    hueseq  = c( seq(5,25,by=2.5), seq(30,75,by=5) )
    
    for( hue in hueseq )
        {
        mess = sprintf( "-------------------  Hue %s  ---------------------\n", HueStringFromNumber(hue) )
        cat( mess )
        
        HVC = cbind( hue, runif(count,min=2,max=8), runif(count,min=1,max=8) )
        
        dat = roundHVC( HVC, book='soil' )
        
        #   hopefully all the rounded hues are hue
        HVCrnd  = HVCfromMunsellName( dat$MunsellRnd )
        maskbad = HVCrnd[ ,1] != hue
        
        if( any(maskbad) )
            {
            print( dat[maskbad, ] )
            countbad    = countbad + sum(maskbad)
            
            # return(FALSE)
            }
        }
        
    #   calculate fraction of bad ones
    fractionbad = countbad / ( count*length(hueseq) )
        
        
    cat( "bad count: ", countbad, " of ",  count*length(hueseq),  "    Fraction bad = ", fractionbad, '\n' )
    
    return( fractionbad < 0.02 )
    }
    
    
if( ! testRounding() )    stop( "testRounding() failed !\n" )


cat(  "Passed rounding test !\n" )
