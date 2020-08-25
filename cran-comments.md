## Test environments
* local R installation, R 4.0.0
* ubuntu 16.04 (on travis-ci), R 4.0.0
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Rhub results 

Aside from the testing CI environments, I also submitted and checked the package compiled and tested 
successfully on the following Rhub platforms:

* debian-clang-devel
* debian-gcc-release
* macos-highsierra-release-cran
* solaris-x86-patched
* ubuntu-gcc-release
* windows-x86_64-release

The solaris build recorded a NOTE about a compiler flag that seemed unrelated to my package. While I do have a Makevars
and Makevars.win to ensure the C++11 standard is used and the headers are included appropriately, I do not make any changes
to the default compiler flags. 