# NEWS for **munsellinterpol** package

### Changes for version 2.5-1  [2019-07-20]

* in testing, allow limit of 5 failures in `testOptimals()`, which started to happen on **fedora** platform (in addition to **solaris**)


### Changes for version 2.4-1  [2019-05-16]

* restored missing .R files in folder /inst/doc


### Changes for version 2.3-1  [2019-04-16]

* for the function `IsWithinMacAdamLimits()`, switched to a zonohedral representation of the color solids
* moved `rootSolve` from dependency to import
* removed function `xyz2srgb()`
* removed dependency `geometry`
* removed folder `demo`


### Changes for version 2.2-1  [2019-01-27]

* added fix to track a change in package 'spacesRGB'
* added Color of the Year for 2019 to User Guide
* moved most startup code from .onAttach() to .onLoad()


### Changes for version 2.1-3  [2018-07-22]

* added the Nickerson Color Difference formula
* refactored by putting XYZ-related functions in new package 'spacesXYZ'
* added precomputed CATs for Illuminant D65  <->  Illuminant C  (for sRGB <-> Munsell conversion)
* fixed issue with no long double (noLD) test


### Changes for version 2.0.1  [2018-06-06]

* the 2 core functions: MunsellToxyY()  and  xyYtoMunsell()  have been re-written, and have new options
* added function `ColorBlockFromMunsell()` for ISCC-NBS color names and centroids
* added functions `plotLociHC()` and `plotPatchesH()`
* added support for wide-gamut Adobe RGB
* added support for user-defined RGB spaces
* improved set of optimal xyY's for Illuminant C
* new munsellinterpol User Guide
* 2 new vignettes
* Imports package `rootSolve`, for the inverse mapping xyY -> HVC
* added extensive test suite
