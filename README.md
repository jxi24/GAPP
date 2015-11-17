#This is a README file.
# Kirtimaan: Attempting to change code to c++
#1. Make code compileable without gfortran
#2. Reorganize directory structure and makefile
#    make file  is now inside bin/makefile

##################################################
############Dependencies #########################
##################################################
Currently using 
1. gcc version 5.1.1 20150618 (Red Hat 5.1.1-4) (GCC)
2. BOOST LIBRARIES
3. ROOT 5.34/32
###################################################
# Description of the structure of the code ########
###################################################
All source files are in directory src/
Fortran codes are in src/F/ and src/F/mdls/
Following is the description of the code also avaiable at
http://arxiv.org/pdf/hep-ph/0005084v1.pdf
###################################################
## Fortran code
###################################################
##Main programs from original GAPP code:
1. fit.f		: contains a simple call to minuit for SM fit
2. mh.f			: Computes proabability distribution function for Mhiggs:(Relevant for pre-higgs discovery era)
3. plotter.f		: Plotting program
4. zplotter.f		: Plotting program
5. zplotter2.f		: Plotting program
####################################################
##Additional main programs from Kai Schmidt
6. average.f 		: Not sure: seems to be independent code that calculates some sort of average-- need to check if it is being called by the scripts
7. sort.f		: Not sure: seems to be independent code that sorts data in a particular file-- need to check if it is being called by the scripts
8. sinth_zprime.f	:
9. mhplot.f		:
10. 221mt.f		:
11. 221mh.f		:
12. 221plotter.f	:
13. plotter2.f		:
14. plotter3.f		:
15. plot_ma.f		:
#####################################################
##Remaing code
13. chi2.f		: contains subroutine 1. fcn which:
						a. Defines constants and flags
						b. Initializes parts of the 1-loop package FF
						c. Makes final call to subroutine values in main.f
					and a function chi2.f which:
						a. Calculates the chi2. The user can make changes to the data , values, errors and correlation coeffs here.
14. main.f 		: Contains subroutine values that 
15. bsgamm.f		: Calculates b to s gamma
16. sin2th.f		: Calculates \sin^{2}\theta_{W}
17. alfahat.f 		: Calculates RGE of \alpha to order (\alpha \alpha_{s}^{3})
18. lep100.f		: Calculates Z-pole observables, decay widths etc.
19. bvertex.f 		: Calculates the top mass corrections to Z-> bb decay
20. masses.f		: Running masses of quarks to 3 loop order
21. alfas.f		: Strong coupling constant at 4 loop
22.

						



