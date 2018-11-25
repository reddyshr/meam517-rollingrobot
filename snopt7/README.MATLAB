================================
SNOPT mex Files in $SNOPT/matlab
================================

The distribution includes several precompiled mex files for SNOPT and
its QP solver SQOPT

  % $SNOPT/matlab/precompiled/snoptmex.mexa64
  % $SNOPT/matlab/precompiled/snoptmex.mexmaci64
  % $SNOPT/matlab/precompiled/snoptmex.mexw64
  % $SNOPT/matlab/precompiled/snoptmex.mexw32

  % $SNOPT/matlab/precompiled/sqoptmex.mexa64
  % $SNOPT/matlab/precompiled/sqoptmex.mexmaci64
  % $SNOPT/matlab/precompiled/sqoptmex.mexw64
  % $SNOPT/matlab/precompiled/sqoptmex.mexw32


You may also compile your own mex files.  To do so, you provide the
location of MATLAB on your machine using:

  >>  ./configure --with-matlab=/path/to/matlab
  >>  make matlab

By default, the mex file will be compiled using the BLAS routines provided
by SNOPT.  To compile against Matlab's Intel BLAS library, SNOPT should be
configured with 64-bit integers and with the blas flag set to "matlab":

  >>  ./configure --with-matlab=/path/to/matlab --with-blas=matlab --with-64
  >>  make matlab

The mex file will be compiled during the make process and placed in the
$SNOPT/matlab directory.


==========================================
Using the SNOPT mex Files in $SNOPT/matlab
==========================================

 Run Matlab in the directory $SNOPT/matlab.

 Typing
  >> more on
  >> help Contents
from Matlab provides an overview of the package.

 At the Matlab prompt, type

  >> runNPExamples    % to test the NLP (snopt) interface
  >> runQPExamples    % to test the QP  (sqopt) interface

This script sets an appropriate matlab path.

The subdirectory  $SNOPT/matlab/examples  contains sample problems that
demonstrate how to use the various snOpt Matlab interfaces.
Read the "Contents" for more information.

Problems may be run individually from $SNOPT/matlab as follows.

    >> setpath    % (if you haven't already called runAllExamples)

    >> addpath examples
    >> addpath examples/snmain
    >> snoptmain  % (this is one of the examples)
