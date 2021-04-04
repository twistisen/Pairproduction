# Pairproduction
A C++ program to calculate the electron/positron pair production using the Semi-classical approach of Baier et al.
This code was used to produce the results shown in Phys. Rev. D 101, 076017 and any use of this code must be accompanied with a reference to this paper.
The output is the file spectrum.txt where the first column is the positron energy and the proceeding 8 columns the different combinations of spin and relative polarization between laser and incoming high energy photon.
Parameters which you may which to change are in constants.h including number of CPU threads to use in the computation (currently set to 7)
