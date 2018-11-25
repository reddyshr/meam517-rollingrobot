README.AMPL.                                   Elizabeth Wong   26 February 2014


AMPL is a powerful programming language for linear and nonlinear programming.


To enable AMPL support in SNOPT, configure SNOPT with the '--with-ampl' switch.

  >> ./configure --with-ampl
  >> make ampl

This will compile the AMPL source files included in $SNOPT and build an AMPL
library, along with SNOPT.

If you already have the AMPL library (amplsolver.a) on your machine, then you
 give the location in the switch, e.g., if amplsolver.a is in /of/lib/, then

  >> ./configure --with-ampl=/location/of/lib

The configure script will search for amplsolver.a in the directory
/location/of/lib/.
