## Test environments
* ubuntu 14.04 (on travis-ci), R 3.4.0
* local windows 10 machine, R 3.4.3
* win-builder

## R CMD check results
There were no ERRORs or WARNINGs. 

When run locally on windows, there was 1 NOTE:

File 'circglmbayes/libs/x64/circglmbayes.dll':
  Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'

It is good practice to register native routines and to disable symbol
search.


This is a known false positive on Windows, as registration is in fact taken care of in RcppExports.R and this NOTE is not happening on ubuntu.


