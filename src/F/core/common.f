      dimension q(0:9),i3(0:9),nc(0:9),v(0:9),a(0:9),v2(0:9),a2(0:9),
     .          eps2_L(0:9),eps2_R(0:9),
     .          asgrid(0:2000),an(4,0:2000),bn(4,0:2000)

C Kai003 ----
C     Dimensions of the additional couplings.
C     i3r = T_R^3, charge under new SU(2).
C     u1xr = U_R, charge of right-handed reps. under U_X(1). 
C     u1xl = U_L, charge of left-handed reps. under U_X(1).
C -----------
      dimension i3r(0:9),u1xr(0:9),u1xl(0:9)
C -----------

C Kai020 ----
C     Switches that allow to print only the pulls of the fitted 
C     observables.
C -----------
      dimension        prtpll(0:80)
C -----------

      double precision alpha,gf,alfas0,asgrid,mu0,dahad3,pol4hf,an,bn,
     .                 pi1,pi2,zeta2,zeta3,zeta4,zeta5,zeta6,zeta7
      double precision q,i3,nc,v,a,alphat,sinhat,coshat,rhonc

C Kai004 ----
C     Additional couplings.
C -----------
      double precision i3r,u1xr,u1xl
C -----------

      double precision mz,mh,mt,mb,mc,ms,md,mu,mtau,mmu,me,mw
      double precision mw2,mz2,mt2,mh2
      double precision ratzw2,rathw2,rathz2,rattw2,rattz2,ratth2
      double precision SWpar,Spar,Tpar,Upar,Brho,Bkappa,Rpar,Zpar,
     .                 delrwh,delrzh,delr,rhohat,rho2
      double precision mzp,mzp2,mz02,sinth,sinth2,costh,costh2,ratg21,
     .                 rhoeff,rhoezp,rhozzp,v2,a2,eps2_L,eps2_R,ratgRL
      double precision prob,mtp,sigma

      logical          flagmr,ffermi,fa2mt4,fa2mt2,fa2mt0,fla2im,
     .                 fasmt2,fasmt0,fas2mt,flagmf,fpolew,flgech,fobliq,
     .                 f4lqcd,falas2,fbayes,flagmh,flagmt,flagmc,flagS,
     .                 flagT,flgrho,fkappa,fzprim,fsinth,fwrite,flprob,
     .                 fhiggs,flagal,fsplot

C Kai025 ----
C     Flags for NP parameters and 'fit221' output
C -----------
      logical          flgfitx, flgtph, flgs2b, flfout
C -----------



      integer          flgblm,flagzp,zprime

C Kai001 ----
C     Additional parameters in the G(221) models. 
C     fitx = u^2 / v^2
C     fitsph = sin^2(phi)
C     fitcph = cos^2(phi)
C     fits2b = sin^2(2beta)
C -----------
      double precision fitx, fitsph, fitcph, fittph,fits2b,
     .                 kkss, kkcc, kkem
C -----------

C Kai020 ----
C     Switch to be able to change between different model types.
C -----------
      integer          modtype
C -----------

C Kai022 ----
C     Triggers to get the initial pulls of the fitted observables.
C -----------
      integer          plltr, prtpll
C -----------

C Kai027 ----
C     Best fit values for internal use in the Fortran files.
C -----------
      double precision bestlnx, besttph, bests2b
C -----------


      common /inputs/  alpha,gf,alfas0,asgrid,mu0,dahad3,pol4hf,an,bn,
     .                 pi1,pi2,zeta2,zeta3,zeta4,zeta5,zeta6,zeta7
C Kai005 ----
C     Update of the common block 'coupls'
C -----------
      common /coupls/  q,i3,nc,v,a,alphat,sinhat,coshat,rhonc,
     .                 i3r,u1xr,u1xl
C -----------
      common /masses/  mz,mh,mt,mb,mc,ms,md,mu,mtau,mmu,me,mw
      common /mass2/   mw2,mz2,mt2,mh2
      common /ratios/  ratzw2,rathw2,rathz2,rattw2,rattz2,ratth2
      common /obliqe/  SWpar,Spar,Tpar,Upar,Brho,Bkappa,Rpar,Zpar,
     .                 delrwh,delrzh,delr,rhohat,rho2
      common /zprime/  mzp,mzp2,mz02,sinth,sinth2,costh,costh2,ratg21,
     .                 rhoeff,rhoezp,rhozzp,v2,a2,eps2_L,eps2_R,ratgRL
      common /limits/  prob,mtp,sigma,zprime
      common /flags/   flagmr,flgblm,ffermi,fa2mt4,fa2mt2,fa2mt0,fla2im,
     .                 fasmt2,fasmt0,fas2mt,flagmf,fpolew,flgech,fobliq,
     .                 f4lqcd,falas2,fbayes,flagzp,flagmh,flagmt,flagmc,
     .                 flgrho,fkappa,flagS,flagT,fzprim,fsinth,fwrite,
     .                 flprob,fhiggs,flagal,fsplot

C Kai002 ----
C     Common block for the new G(221) parameters.
C -----------
      common /fit221/ fitx, fitsph, fitcph, fittph, fits2b,
     .                kkss, kkcc, kkem
C -----------

C Kai020 ----
C     Common block for switches and triggers.
C -----------
      common /fitswtr/ modtype, plltr, prtpll
C -----------

C Kai025 ----
C     Common block for new flags.
C -----------
      common /fitflg/  flgfitx, flgtph, flgs2b, flfout
C -----------

C Kai027 ----
C     Common block for internal best fit values.
C -----------
      common /fitbst/  bestlnx, besttph, bests2b
C -----------