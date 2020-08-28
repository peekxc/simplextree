## Resubmission
This is a resubmission. In this version I have:

* Fixed the Authors@R field again so that the maintainer matches the author field. 

* Fixed the Authors@R field to be a vector of person objects as opposed to as.person

* Added Howard Hinnant to the author list in the DESCRIPTION file as a copyright holder, and fixed the LICENSE.

* Added the doi linking the paper this package was inspired by in the Description tag of the DESCRIPTION file.

* Removed unnecessary dontrun sections in the examples.

* Removed roxygen-generated documentation pages from unexported functions. 

## Test environments
* Local R installation (64-bit x86_64-apple-darwin17.0), R 4.0.0
* Ubuntu 16.04 (on travis-ci), R 4.0.2
* MacOS 10.13 (on travis-ci), R 4.0.2 
* Windows Server 2012 R2 x64 (on appveyor), R 4.0.2 

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

The solaris build recorded a NOTE about a compiler flag that is unrelated to my package. While I do have a Makevars
and Makevars.win to ensure the C++11 standard is used and the headers are included appropriately, I do not make any changes
to the default compiler flags. 