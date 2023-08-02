# HYDROdynamic Freeze-Out IntegraL - hydro-foil
Compute integrals on the freeze-out hypersurface. By default, the integration happens parallelly on NTHREADS = number_of_threads_available - 2. This number can be modified prior to compiling in "integrals.cpp". The use of openmp can be completely switched off commenting the OPENMP flag in the makefile

 The typical execution from command line is:
./foil PATH_TO/beta.dat PATH_TO/outputfile

The main file should be modified changing the function to use depending on the calculation.
