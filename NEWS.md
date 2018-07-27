# NEWS for **munsellinterpol** package

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
