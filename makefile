#FLIBS=   /opt/2002/lib/libpacklib.a \
#         /opt/2002/lib/libmathlib.a 
# CJY
FLIBS = -Lcernlib -L/usr/lib64/cernlib/2006-g77/lib -lpacklib -lnsl -lkernlib
#FLIBS = -Lcernlib -L/scratch/kai/masters/cernlib/cernlib/2006b/i686-slc5-gcc34-opt/lib \
        -lpacklib -lkernlib -lpacklib -lnsl
#CJY
COMPILE.f=f77 -w -c -C -g
LINK.f=xlf -g
FFLAGS= -g

fitobj  = fit.o chi2.o alfahat.o sin2th.o lowenergy.o alfas.o masses.o \
          kappa.o rho.o bvertex.o lep100.o pnc.o dis.o nue.o bsgamma.o \
          polylogs.o spence.o ffcli2.o ffxli2.o ffinit.o main.o zprime.o \
          wwidth.o lep200.o taulife.o amu.o sumrules.o sin2th_plot.o

fitobj2  = fit_v2.o chi2_v2.o alfahat.o sin2th.o lowenergy.o alfas.o masses.o \
          kappa.o rho.o bvertex.o lep100.o pnc.o dis.o nue.o bsgamma.o \
          polylogs.o spence.o ffcli2.o ffxli2.o ffinit.o main.o zprime.o \
          wwidth.o lep200.o taulife.o amu.o sumrules.o sin2th_plot.o

mhobjs  = mh.o chi2.o alfahat.o sin2th.o lowenergy.o alfas.o masses.o \
          kappa.o rho.o bvertex.o lep100.o pnc.o dis.o nue.o bsgamma.o \
          polylogs.o spence.o ffcli2.o ffxli2.o ffinit.o main.o zprime.o \
          wwidth.o lep200.o taulife.o amu.o sumrules.o sin2th_plot.o

mhpobj  = mhplot.o chi2.o alfahat.o sin2th.o lowenergy.o alfas.o masses.o \
          kappa.o rho.o bvertex.o lep100.o pnc.o dis.o nue.o bsgamma.o \
          polylogs.o spence.o ffcli2.o ffxli2.o ffinit.o main.o zprime.o \
          wwidth.o lep200.o taulife.o amu.o sumrules.o sin2th_plot.o

sthobj  = chi2.o alfahat.o sin2th.o lowenergy.o alfas.o masses.o \
          kappa.o rho.o bvertex.o lep100.o pnc.o dis.o nue.o bsgamma.o \
          polylogs.o spence.o ffcli2.o ffxli2.o ffinit.o main.o zprime.o \
          wwidth.o sinth_zprime.o lep200.o taulife.o amu.o sumrules.o \
          sin2th_plot.o

zpobjs  = zplotter.o chi2.o alfahat.o sin2th.o lowenergy.o alfas.o masses.o \
          kappa.o rho.o bvertex.o lep100.o pnc.o dis.o nue.o bsgamma.o \
          polylogs.o spence.o ffcli2.o ffxli2.o ffinit.o main.o zprime.o \
          wwidth.o lep200.o taulife.o amu.o sumrules.o sin2th_plot.o

zpobjs2 = zplotter2.o chi2.o alfahat.o sin2th.o lowenergy.o alfas.o masses.o \
          kappa.o rho.o bvertex.o lep100.o pnc.o dis.o nue.o bsgamma.o \
          polylogs.o spence.o ffcli2.o ffxli2.o ffinit.o main.o zprime.o \
          wwidth.o lep200.o taulife.o amu.o sumrules.o sin2th_plot.o

maobjs  = plot_ma.o chi2.o alfahat.o sin2th.o lowenergy.o alfas.o masses.o \
          kappa.o rho.o bvertex.o lep100.o pnc.o dis.o nue.o bsgamma.o \
          polylogs.o spence.o ffcli2.o ffxli2.o ffinit.o main.o zprime.o \
          wwidth.o lep200.o taulife.o amu.o sumrules.o sin2th_plot.o

plotob  = plotter.o chi2.o alfahat.o sin2th.o lowenergy.o alfas.o masses.o \
          kappa.o rho.o bvertex.o lep100.o pnc.o dis.o nue.o bsgamma.o \
          polylogs.o spence.o ffcli2.o ffxli2.o ffinit.o main.o zprime.o \
          wwidth.o lep200.o taulife.o amu.o sumrules.o sin2th_plot.o

plotob2 = plotter2.o chi2.o alfahat.o sin2th.o lowenergy.o alfas.o masses.o \
          kappa.o rho.o bvertex.o lep100.o pnc.o dis.o nue.o bsgamma.o \
          polylogs.o spence.o ffcli2.o ffxli2.o ffinit.o main.o zprime.o \
          wwidth.o lep200.o taulife.o amu.o sumrules.o sin2th_plot.o

plotob3 = plotter3.o chi2.o alfahat.o sin2th.o lowenergy.o alfas.o masses.o \
          kappa.o rho.o bvertex.o lep100.o pnc.o dis.o nue.o bsgamma.o \
          polylogs.o spence.o ffcli2.o ffxli2.o ffinit.o main.o zprime.o \
          wwidth.o lep200.o taulife.o amu.o sumrules.o sin2th_plot.o

abobjs  = ab.o chi2.o alfahat.o sin2th.o lowenergy.o alfas.o masses.o \
          kappa.o rho.o bvertex.o lep100.o pnc.o dis.o nue.o bsgamma.o \
          polylogs.o spence.o ffcli2.o ffxli2.o ffinit.o main.o zprime.o \
          wwidth.o lep200.o taulife.o amu.o sumrules.o sin2th_plot.o

# Kai024 ###############################################################
# Object files needed to compile 221plotter.f
########################################################################

221plotob = 221plotter.o chi2.o alfahat.o sin2th.o lowenergy.o alfas.o \
            masses.o kappa.o rho.o bvertex.o lep100.o pnc.o dis.o nue.o \
            bsgamma.o polylogs.o spence.o ffcli2.o ffxli2.o ffinit.o main.o \
            zprime.o wwidth.o lep200.o taulife.o amu.o sumrules.o \
            sin2th_plot.o

# Kai029 ###############################################################
# Object files needed to compile 221mh.f and 221mt.f
########################################################################

221mhob = 221mh.o chi2.o alfahat.o sin2th.o lowenergy.o alfas.o \
          masses.o kappa.o rho.o bvertex.o lep100.o pnc.o dis.o nue.o \
          bsgamma.o polylogs.o spence.o ffcli2.o ffxli2.o ffinit.o main.o \
          zprime.o wwidth.o lep200.o taulife.o amu.o sumrules.o \
          sin2th_plot.o

221mtob = 221mt.o chi2.o alfahat.o sin2th.o lowenergy.o alfas.o \
          masses.o kappa.o rho.o bvertex.o lep100.o pnc.o dis.o nue.o \
          bsgamma.o polylogs.o spence.o ffcli2.o ffxli2.o ffinit.o main.o \
          zprime.o wwidth.o lep200.o taulife.o amu.o sumrules.o \
          sin2th_plot.o

########################################################################


fit:    $(fitobj)
	$(FC) -o $@ $(FFLAGS) $(fitobj) $(FLIBS)

fit2:   $(fitobj2)
	$(FC) -o $@ $(FFLAGS) $(fitobj2) $(FLIBS)

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

ab:     $(abobjs)
	$(FC) -o $@ $(FFLAGS) $(abobjs) $(FLIBS)

clean:
	rm -r *.o *-d.fit *-t.fit

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
