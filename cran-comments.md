## Resubmission
This is a resubmission. In this version I have:
* updated the reference format in DESCRIPTION file. From format like (2012)(<doi:10.1080/09332480.2014.914768>) to (2012, <doi:10.1080/09332480.2014.914768>).
* explained the acronyms in the DESCRIPTION file.
* explained the acronyms in some of the Rd files.
* replaced \dontrun{} with \donttest{} in the Rd files.
* added \value tags to all the Rd files where \value was missing and explained the returned values.
* removed the examples with un-exported functions.

One of the feedback from the previous submissoin says "You are changing the user's par() settings in your examples. Please reset the settings at the end of your examples." But there's no par() usage in the examples. Could it because there's a parameter named "parH0" and some how been mistaken as "par()"?

## Test environments
* local Ubuntu 16.04, R 3.6.2
* local Mac OS 10.10, R. 3.6.2
* Ubuntu 16.04 (on travis-ci), R-devel, R 3.5.3, R 3.6.2
* win-builder, R-devel

## R CMD check results
0 errors ✔ | 0 warnings ✔ | 0 notes ✔


## Downstream dependencies
There are currently no downstream dependencies for this package
