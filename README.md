# HYDROdynamic Freeze-Out IntegraL - hydro-foil
Compute integrals on the freeze-out hypersurface. By default, the integration happens parallelly on NTHREADS = number_of_threads_available - 2. This number can be modified prior to compiling in "integrals.cpp". The use of OpenMP can be completely switched off by commenting the "OPENMP" flag in the makefile

 The typical execution from command line is:
./foil PATH_TO/beta.dat PATH_TO/outputfile

The main file should be modified changing the function to use depending on the calculation. The functions available at the moment are:

* **void polarization_midrapidity(double pT, double phi, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout);** Computes the mean spin vector for particles emitted at midrapidity for the particle "particle" using the freezeout file "freezeout_sup". For the thermal vorticity, eq. 64 of 2304.02276 is used. Refer also to 2103.10917 for details about the formulae.
The output is written in "fileout" as a table:

| pt | phi | denominator | numerator_varpi | numerator_xi |

* **void polarization_midrapidity_linear(double pT, double phi, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout);** Same as the previous function but uses the linear approximation for the vorticity-induced polarization. This function is faster.



* **void polarization_projected(double pT, double phi, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout);** Similar to "polarization_midrapidity_linear" but computes polarization projecting the gradients tangent to the hypersurface at point x before integrating.  

* **void polarization_rapidity(double pT, double phi, double y_rap, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout);** Same as "polarization_midrapidity", but the table now includes the rapidity "y":

| pt | phi | y | denominator | numerator_varpi | numerator_xi |


* **void spectrum_rapidity(double pT, double phi, double y_rap, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout);** Integrates the thermal spectrum of "particle" using the Cooper-Freye prescription.


