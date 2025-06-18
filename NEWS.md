# NEWS for **munsellinterpol** package


### Version 3.2-0  [2025-06-17]

* added new vignette: Soil Colors
* added arguments `value` and `chroma` to `plotPatchesH()`
* added argument `coeffs` to `NickersonColorDifference()`


### Version 3.1-0  [2025-01-19]

* all logging now done with package **logger**, which is Imported
* updates to User Guide, including 3 new Pantone Colors of the Year
* in function `xyYtoMunsell()`, add new arguments `rtol` and `atol`


### Version 3.0-0  [2022-04-08]

* for conversion functions not involving `XYZ`, ensure that an exact neutral converts to an exact neutral
* for conversion functions involving `Lab` and `Luv`, allow reference white to be assigned by name
* consistent handling of `rownames` in conversion function return values
* for functions `LabtoMunsell` and `LuvtoMunsell` changed `t` to `T` (old spelling is still supported for a limited time)
* to functions `RGBtoMunsell()` and `MunsellToRGB()`, added argument `which` 
* in the User Guide, added 3 new Pantone Colors of the Year
* in links to other packages from man pages, added package names
* added hue circle figure to the PDF manual


### Version 2.8-2  [2022-03-03]

* in one of the `.Rd` files, replaced `<center>` tag by the HTML5 equivalent
* removed some (possibly) invalid URLs


### Version 2.7-1  [2021-12-09]

* fixed `bibliography.bib` to be compatible with `pandoc` v. 2.16.2
* fixed some stale URLs in `bibliography.bib` and man pages
* removed superfluous argument in a call to `sprintf()`


### Changes for version 2.6-1  [2020-02-02]

* fixed error in a test regarding `identical()` vs `all.equal()`  (related to noLD)
* suppressed warning from `stats::regularize.values` by adding argument `ties=min` to `stats::splinefun()`


### Changes for version 2.5-1  [2019-07-20]

* in testing, allow limit of 5 failures in `testOptimals()`, which started happening on **fedora** platform (in addition to **solaris**)


### Changes for version 2.4-1  [2019-05-16]

* restored missing .R files in folder /inst/doc


### Changes for version 2.3-1  [2019-04-16]

* for the function `IsWithinMacAdamLimits()`, switched to a zonohedral representation of the color solids
* moved package `rootSolve` from dependency to import
* removed function `xyz2srgb()`
* removed dependency `geometry`
* removed folder `demo`


### Changes for version 2.2-1  [2019-01-27]

* added fix to track a change in package `spacesRGB`
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
