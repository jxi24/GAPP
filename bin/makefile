#FLIBS=   /opt/2002/lib/libpacklib.a \
#         /opt/2002/lib/libmathlib.a 
# CJY
#FLIBS = -Lcernlib -L/usr/lib64/cernlib/2006-g77/lib -lpacklib -lnsl -lkernlib
#FLIBS = -Lcernlib -L/scratch/kai/masters/cernlib/cernlib/2006b/i686-slc5-gcc34-opt/lib \
        -lpacklib -lkernlib -lpacklib -lnsl
#CJY
SOURCE=../src/F/core/:../src/F/driver_code/
BUILD=../build/

VPATH=$(SOURCE)

FC=gfortran
FCC=gfortran -w -c -C -g
LINK.f=xlf -g
FFLAGS= -g


commonobj = $(BUILD)alfahat.o $(BUILD)sin2th.o $(BUILD)lowenergy.o $(BUILD)alfas.o $(BUILD)masses.o \
            $(BUILD)kappa.o $(BUILD)rho.o $(BUILD)bvertex.o $(BUILD)lep100.o $(BUILD)pnc.o $(BUILD)dis.o $(BUILD)nue.o $(BUILD)bsgamma.o \
            $(BUILD)polylogs.o $(BUILD)spence.o $(BUILD)ffcli2.o $(BUILD)ffxli2.o $(BUILD)ffinit.o $(BUILD)main.o $(BUILD)zprime.o \
            $(BUILD)wwidth.o $(BUILD)lep200.o $(BUILD)taulife.o $(BUILD)amu.o $(BUILD)sumrules.o $(BUILD)sin2th_plot.o

fitobj  = $(BUILD)fit.o $(BUILD)chi2.o $(commonobj)

fitobj2  = $(BUILD)fit_v2.o $(BUILD)chi2_v2.o $(commonobj)

mhobjs  = $(BUILD)mh.o $(BUILD)chi2.o $(commonobj)

mhpobj  = $(BUILD)mhplot.o $(BUILD)chi2.o $(commonobj)

sthobj  = $(BUILD)chi2.o $(BUILD)sinth_zprime.o $(commonobj)

zpobjs  = $(BUILD)zplotter.o $(BUILD)chi2.o $(commonobj)

zpobjs2 = $(BUILD)zplotter2.o $(BUILD)chi2.o $(commonobj)

maobjs  = $(BUILD)plot_ma.o $(BUILD)chi2.o $(commonobj)

plotob  = $(BUILD)plotter.o $(BUILD)chi2.o $(commonobj)

plotob2 = $(BUILD)plotter2.o $(BUILD)chi2.o $(commonobj)

plotob3 = $(BUILD)plotter3.o $(BUILD)chi2.o $(commonobj)

abobjs  = $(BUILD)ab.o $(BUILD)chi2.o $(commonobj)

# Kai024 ###############################################################
# Object files needed to compile 221plotter.f
########################################################################

221plotob = $(BUILD)221plotter.o $(BUILD)chi2.o $(commonobj)

# Kai029 ###############################################################
# Object files needed to compile 221mh.f and 221mt.f
########################################################################

221mhob = $(BUILD)221mh.o $(BUILD)chi2.o $(commonobj)

221mtob = $(BUILD)221mt.o $(BUILD)chi2.o $(commonobj)

########################################################################


fit:    $(fitobj)
	$(FC) -o $@ $(FFLAGS) $(fitobj) $(FLIBS)

## Could not find file fit_f2.f ==> Check this 
#fit2:   $(fitobj2)
#	$(FC) -o $@ $(FFLAGS) $(fitobj2) $(FLIBS)

mh:     $(mhobjs)
	$(FC) -o $@ $(FFLAGS) $(mhobjs) $(FLIBS)

mhplot: $(mhpobj)
	$(FC) -o $@ $(FFLAGS) $(mhpobj) $(FLIBS)

sinth:  $(sthobj)
	$(FC) -o $@ $(FFLAGS) $(sthobj) $(FLIBS)

zplot:  $(zpobjs)
	$(FC) -o $@ $(FFLAGS) $(zpobjs) $(FLIBS)

zplot2: $(zpobjs2)
	$(FC) -o $@ $(FFLAGS) $(zpobjs2) $(FLIBS)

maplot: $(maobjs)
	$(FC) -o $@ $(FFLAGS) $(maobjs) $(FLIBS)

plot:   $(plotob)
	$(FC) -o $@ $(FFLAGS) $(plotob) $(FLIBS)

plot2:  $(plotob2)
	$(FC) -o $@ $(FFLAGS) $(plotob2) $(FLIBS)

plot3:  $(plotob3)
	$(FC) -o $@ $(FFLAGS) $(plotob3) $(FLIBS)

#ab:     $(abobjs)
#	$(FC) -o $@ $(FFLAGS) $(abobjs) $(FLIBS)

# Kai024 ###############################################################
# Make command for the program ttoplot in 221plotter.f 
########################################################################

221plot:$(221plotob)
	$(FC) -o $@ $(FFLAGS) $(221plotob) $(FLIBS)

# Kai029 ###############################################################
# Make command for the programs ttomh and ttomt in 221mh.f and 221mt.f
########################################################################

221mh:  $(221mhob)
	$(FC) -o $@ $(FFLAGS) $(221mhob) $(FLIBS)

221mt:  $(221mtob)
	$(FC) -o $@ $(FFLAGS) $(221mtob) $(FLIBS)

########################################################################

########################################################################

$(BUILD)%.o: $(SRC)%.f
	$(FCC) $< -o $@ $(FFLAGS) $(FLIBS)

LIST= fit mh zplot zplot2 maplot plot plot2 plot3  221plot 221mh 221mt

all: $(LIST)

clean:
	rm  -f $(BUILD)*.o 
	rm -f *-d.fit *-t.fit
cleanall:
	rm -f $(LIST)
	rm  -f $(BUILD)*.o 
	rm -f *-d.fit *-t.fit


