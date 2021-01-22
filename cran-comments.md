## Test environments
* Ubuntu 16.04.6 LTS (on travis-ci), R 4.0.2
* local windows 10 machine, 4.0.3
* win-builder (release & devel)

## R CMD check results
There were no ERRORs or WARNINGs. 

When run locally on windows, there was 1 NOTE:

    > checking compiled code ... NOTE
      Note: information on .o files for i386 is not available
      Note: information on .o files for x64 is not available
      File 'C:/Dropbox/Research/circglmbayes.Rcheck/circglmbayes/libs/i386/circglmbayes.dll':
        Found 'abort', possibly from 'abort' (C), 'runtime' (Fortran)
        Found 'exit', possibly from 'exit' (C), 'stop' (Fortran)
        Found 'printf', possibly from 'printf' (C)

It seems like this is a known false positive (e.g. https://stackoverflow.com/questions/64402688/information-on-o-files-for-x64-is-not-available-note-on-r-package-checks-using), and seems not to be a result of what is done by this package. 
