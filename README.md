# ASP

Active-set procedures for Basis Pursuit and related  problems.

## Authors

- Michael Friedlander <michael.friedlander@ubc.ca>
- Michael Saunders <saunders@stanford.edu>


## Install and test BPdual

    cd <asp directory>   % change to the asp directory
    addpath(pwd)         % add it to the path
    as_setup             % compile mex interfaces
    cd tests
    demo_bpdn            % run a demo

## Reference

Michael Friedlander and Michael Saunders. [A dual active-set
quadratic programming method for finding sparse least-squares
solutions](https://friedlander.io/files/pdf/bpprimal.pdf), DRAFT Technical Report, Dept of Computer Science,
University of British Columbia, July 30, 2012; updated April 13, 2019.
