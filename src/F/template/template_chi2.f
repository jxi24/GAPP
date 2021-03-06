      subroutine fcn(npar,grad,fval,xval,iflag,chi2)

      implicit none
      integer npar,iflag,i,j
      double precision grad(npar),fval,xval(npar),chi2,det,dgamma
      double precision smval(80),pull(80),emat(6,6),dummy
      integer ii2

      include 'common.f'
      external chi2

      if (iflag.eq.1) then
         alpha = 1/137.03599911d0
            gf =     1.16637d-5

         ms   =  120.000000d-3  ! +21 -26 MeV (179) in hep-ph/0507078
         md   =    5.440000d-3  ! light quark masses at mu = mtau 
         mu   =    2.180000d-3
         mtau =    1.776990d0
         mmu  =  105.658357d-3
         me   =  510.998902d-6

         q(0) =  0.d0
         q(1) = -1.d0
         q(2) = -1.d0
         q(3) = -1.d0
         q(4) =  2.d0/3
         q(5) =  2.d0/3
         q(6) =  2.d0/3
         q(7) = -1.d0/3
         q(8) = -1.d0/3
         q(9) = -1.d0/3
         
         i3(0) =  1.d0/2
         i3(1) = -1.d0/2
         i3(2) = -1.d0/2
         i3(3) = -1.d0/2
         i3(4) =  1.d0/2
         i3(5) =  1.d0/2
         i3(6) =  1.d0/2
         i3(7) = -1.d0/2
         i3(8) = -1.d0/2
         i3(9) = -1.d0/2
         
         nc(0) = 3.d0
         nc(1) = 1.d0
         nc(2) = 1.d0
         nc(3) = 1.d0
         nc(4) = 3.d0
         nc(5) = 3.d0
         nc(6) = 3.d0
         nc(7) = 3.d0
         nc(8) = 3.d0
         nc(9) = 3.d0

C Kai021 ---------------------------------
C  Model type switch.
        
         modtype = 0
C ----------------------------------------

C Kai007 ---------------------------------
C  Couplings fetched from external file.
        
         include 'mdls/MODELNAME.f'     
C ----------------------------------------

C Kai022 ---------------------------------
C  Pull trigger.
        
         plltr = 1
C ----------------------------------------

         
            pi1 = 4.d0*datan(1.d0)
            pi2 = pi1**2
          zeta2 = pi2/6
          zeta3 = 1.2020569031595942853997382d0
          zeta4 = pi2**2/90
          zeta5 = 1.0369277551433699263313655d0
          zeta6 = pi2**3/945 
          zeta7 = 1.0083492773819228268397975d0
         pol4hf = 0.5174790616738993863307582d0

**********************************************************************
*                                                                    *
*                  FLAGS: (default: .true.)                          *
*                  ------                                            *
*                                                                    *
*   Marciano-Rosner decoupling?                                      *
*                                                                    *
         flagmr = .true.
*                                                                    *
*                                                                    *
*   top pole mass: 0 = 3-loop, 1 = 3-loop "resummed", 2 = 2-loop BLM *
*                                                                    *
         flgblm = 0
*                                                                    *
*                                                                    *
*   Effective charge optimization for O(alpha_s**4) term in Z-width? *
*                                                                    *
         flgech = .true.
*                                                                    *
*                                                                    *
*   4-loop running of alpha_s                                        *
*                                                                    *
         f4lqcd = .true.
*                                                                    *
*                                                                    *
*   Use of Fermi constant in O(GF mt**2) corrections?                *
*                                                                    *
         ffermi = .true.
*                                                                    *
*                                                                    *
*   2-loop O(alpha**2 mt**4) corrections?                            *
*                                                                    *
         fa2mt4 = .true.
*                                                                    *
*                                                                    *
*   2-loop O(alpha**2 mt**2) corrections?                            *
*                                                                    *
         fa2mt2 = .true.
*                                                                    *
*                                                                    *
*   2-loop O(alpha**2 mt**0) corrections (photonic only)?            *
*                                                                    *
         fa2mt0 = .true.
*                                                                    *
*                                                                    *
*   2-loop O(alpha**2) improvement (imaginary parts)?                *
*                                                                    *
         fla2im = .true.
*                                                                    *
*                                                                    *
*   2-loop O(alpha alpha_s mt**2) corrections?                       *
*                                                                    *
         fasmt2 = .true.
*                                                                    *
*                                                                    *
*   2-loop O(alpha alpha_s mt**0) corrections?                       *
*                                                                    *
         fasmt0 = .true.
*                                                                    *
*                                                                    *
*   3-loop O(alpha alpha_s**2 mt**2) corrections?                    *
*                                                                    *
         fas2mt = .true.
*                                                                    *
*                                                                    *
*   3-loop O(alpha alpha_s**2 mt**0) corrections?                    *
*                                                                    *
         falas2 = .true.
*                                                                    *
*                                                                    *
*   light fermion mass effects?                                      *
*                                                                    *
         flagmf = .true.
*                                                                    *
*                                                                    *
*   Higgs direct search results?                                     *
*                                                                    *
         fhiggs = .false.
*                                                                    *
*                                                                    *
*   Bayesian improvement of M_H limit?                               *
*                                                                    *
         fbayes = .false.
*                                                                    *
*                                                                    *
*   oblique parameters for new physics?                              *
*                                                                    *
         fobliq = .true.
*                                                                    *
*                                                                    *
*   no extra Z boson:                                    flagzp = 0  *
*   Z_chi defined by SO(10) -->  SU(5) X U(1)_chi:       flagzp = 1  *
*   Z_psi defined by    E_6 --> SO(10) X U(1)_psi:       flagzp = 2  *
*   Z_eta defined as sqrt(3/8) Z_chi - sqrt(5/8) Z_psi:  flagzp = 3  *
*   Z_LR: SU(2)_R X U(1)_(B-L) --> U(1)_Y X U(1)_LR:     flagzp = 4  *
*   sequential Z^prime:                                  flagzp = 5  *
*   Z^string:                                            flagzp = 6  *
*   the general E_6 boson:                               flagzp = 7  *
*   the general family-universal Z^prime boson:          flagzp = 8  *
*   the model independent Z^prime boson:                 flagzp = 9  *
*   model independent Z^prime in v,a parametrization:    flagzp = 10 *
*   the general E_6 boson in alpha,beta parametrization: flagzp = 11 *
*   the twisted E_6 boson:                               flagzp = 12 *
*   seesaw neutrino mass related (Adhikari, Erler, Ma):  flagzp = 13 *
*                                                                    *
         flagzp = 0
*                                                                    *
*                                                                    *
**********************************************************************

C  Input of functions A_n from grid.out for tau lifetime and g-2:

         if (f4lqcd.eqv..true.) then
            open (3,file='../src/F/dat/grid4loop.dat',status='old')
            open (4,file='../src/F/dat/grid4amu.dat' ,status='old')
         else
            open (3,file='../src/F/dat/grid3loop.dat',status='old')
            open (4,file='../src/F/dat/grid3amu.dat' ,status='old')
         endif

         do 5 i = 0,2000
            read(3,10) asgrid(i),(an(j,i),j=1,4) 
            read(4,10) asgrid(i),(bn(j,i),j=1,4) 
 5       continue

         close(3)
         close(4)
         
 10      format (f6.4,4(1x,f17.14))
         
         call ffinit
         print*,"====================================="
      endif
      fval = chi2(xval,npar,smval,pull)
      prob = fval
c      write(*,*) ("xval[",ii2,"]= ",xval(ii2), ii2=1,npar )
c      print*,"------------" 
c      print*,npar,fval

***********IF Minuit finishes minimizing then do the followin **********
      if (iflag.eq.3) then 
C Kai026 ---------------------------------
C  Store the best fit values of the NP parameters. 
C ----------------------------------------

       open (1,
     $ file='../221plots/bestfit_MODELNAME.tmp'
     $ ,status='unknown')

      write(1,'(F10.5)') xval(2)
      write(1,'(F10.5)') xval(7)
      write(1,'(F10.5)') xval(28)
      write(1,'(F10.5)') xval(29)
      write(1,'(F10.5)') xval(30)

      close(1)

      bestlnx = xval(28)
      besttph = xval(29)
      bests2b = xval(30)

C ----------------------------------------

         mtp = smval(31)
         if (fwrite.eqv..true.) then
            write(7,*) 'eps2_L(e)   = ', eps2_L(1)
            write(7,*) 'eps2_R(e)   = ', eps2_R(1)
            write(7,*) 'eps2_L(u,d) = ', eps2_L(4)
            write(7,*) 'eps2_R(u)   = ', eps2_R(4)
            write(7,*) 'eps2_R(d)   = ', eps2_R(7)
            write(7,*) 'eps2_L(tau) = ', eps2_L(3)
            write(7,*) 'eps2_R(tau) = ', eps2_R(3)
            write(7,*) 'eps2_L(t,b) = ', eps2_L(9)
            write(7,*) 'eps2_R(b)   = ', eps2_R(9)
            write(7,*)
            
            write(7,*) 'v2(nu)  = ', v2(0)
            write(7,*) 'a2(nu)  = ', a2(0)
            write(7,*) 'v2(e)   = ', v2(1)
            write(7,*) 'a2(e)   = ', a2(1)
            write(7,*) 'v2(u)   = ', v2(4)
            write(7,*) 'a2(u)   = ', a2(4)
            write(7,*) 'v2(d)   = ', v2(7)
            write(7,*) 'a2(d)   = ', a2(7)
            write(7,*) 'v2(tau) = ', v2(3)
            write(7,*) 'a2(tau) = ', a2(3)
            write(7,*) 'v2(b)   = ', v2(9)
            write(7,*) 'a2(b)   = ', a2(9)
            write(7,*)

            write(7,100) (i,smval(i),pull(i),i=1,80)
            if (fsplot.eqv..true.) call s2plot
            call values
         endif

C Kai022 ---------------------------------
C  Initial pulls
C ----------------------------------------

      if (flfout.eqv..true.) then  
         do 900 i = 1, 80
            if (prtpll(i).eq.1) then
               write(50,FMT='(I3,A,F8.4)') 100 + i,'\t', pull(i)
            endif
 900     continue
        close(50)
C ----------------------------------------

C Kai024 ---------------------------------
C  Write result to summary.out
C ----------------------------------------
        
        write(51,FMT='(A)') '// MODELNAME ////////////////////////'
     
C       close(51)

      endif
C ----------------------------------------

C---- remove this part of code when using with C++ throught the preprocessor defintion 
#ifndef CPP_FLAG
         if (flprob.eqv..true.) prob = dexp( - fval/2)
         if (fbayes.eqv..true.) then
            call mnemat(emat,6)
            det = emat(1,1)*(emat(2,2)*emat(3,3) - emat(2,3)*emat(3,2))
     .          - emat(1,2)*(emat(2,1)*emat(3,3) - emat(2,3)*emat(3,1))
     .          + emat(1,3)*(emat(2,1)*emat(3,2) - emat(2,2)*emat(3,1))
            prob = prob*dsqrt(det)
         endif
#endif         
      endif
      
 100  format(' observable ',i2,'): ',f12.7,f10.3)

      return
      end
      
************************************************************************

      double precision function chi2(xval,npar,smval,pull)
   
      implicit none
      integer i,j,npar,ima,ssteps
      double precision gamma(0:11),sigmah,R(0:9),sin2te(0:9),alr(0:10)
      double precision afb(0:10),dkapse,C1u,C1d,C2u,C2d,sgcenu,therr2
      double precision xval(npar),smval(80),value(80),error(80),pull(80)
      double precision epsu_L,epsd_L,epsu_R,epsd_R,dummy,fac,lambdg,amu
      double precision ccls(9,9),cchf(15:20,15:20),ccdis(49:55,49:55)
      double precision ccslc(21:24,21:24),chiggs,alphas,xs200(2,0:10)
      double precision afb200(2,0:10),errorc(0:8),errorb(0:8),ccmbc
      double precision cctau,ccmgw,ccntv,ccnue,cchad,ccmt(31:41,31:41)
      double precision momntc(8),momntb(8),sumrlc(0:8),sumrlb(0:8)
      double precision gammaw(7),shat
      complex*16 kappa(0:9),rho(0:9)
      logical flcorr,flgzph
      common   /ma/ ima,ssteps
      integer ii2
      include 'common.f'

***********************************************************************
*                                                                     *
*                   list of observables:                              *
*                   --------------------                              *
*                                                                     *
*       Z lineshape:                                                  *
*                                                                     *
*    1) M_Z                                                           *
*    2) Gamma_Z                                                       *
*    3) sigma_had                                                     *
*    4) R_e                                                           *
*    5) R_mu                                                          *
*    6) R_tau                                                         *
*    7) A^FB (e)                                                      *
*    8) A^FB (mu)                                                     *
*    9) A^FB (tau)                                                    *
*                                                                     *
*       other LEP 1 measurements:                                     *
*                                                                     *
*   10) P (tau)                                                       *
*   11) P^FB (tau)                                                    *
*   12) sin^2 theta^eff_e (Q_FB)                                      *
*   13) A^FB (s) (DELPHI + OPAL)                                      *
*   14) R_d,s/(R_d + R_u + R_s)                                       *
*                                                                     *
*       LEP 1 and SLC heavy flavor:                                   *
*                                                                     *
*   15) R_b                                                           *
*   16) R_c                                                           *
*   17) A^FB (b)                                                      *
*   18) A^FB (c)                                                      *
*   19) A_LR^FB (b)                                                   *
*   20) A_LR^FB (c)                                                   *
*                                                                     *
*       other SLD asymmetries:                                        *
*                                                                     *
*   21) A_LR (hadrons)                                                *
*   22) A_LR (leptons)                                                *
*   23) A_LR^FB (mu)                                                  *
*   24) A_LR^FB (tau)                                                 *
*   25) A_e (Q_LR)                                                    *
*   26) A_LR^FB (s)                                                   *
*                                                                     *
*       W properties:                                                 *
*                                                                     *
*   27) M_W     (LEP)                                                 *
*   28) Gamma_W (LEP)                                                 *
*   29) M_W     (Tevatron)                                            *
*   30) Gamma_W (Tevatron)                                            *
*                                                                     *
*       top quark mass:                                               *
*                                                                     *
*   31) m_t (pole) lepton + jets (CDF I)                              *
*   32) m_t (pole) dilepton      (CDF I)                              *
*   33) m_t (pole) all hadronic  (CDF I)                              *
*   34) m_t (pole) lepton + jets (D0 I)                               *
*   35) m_t (pole) dilepton      (D0 I)                               *
*   36) m_t (pole) lepton + jets (CDF II)                             *
*   37) m_t (pole) dilepton      (CDF II)                             *
*   38) m_t (pole) all hadronic  (CDF II)                             *
*   39) m_t (pole) l+j JES free  (CDF II)                             *
*   40) m_t (pole) lepton + jets (D0 II)                              * 
*   41) m_t (pole) dilepton      (D0 II)                              *
*   42) m_t (pole) all hadronic  (D0 I)   ! hep-ex/0410086            *
*   42) m_t (pole) TOTAL TEVATRON average incl. 0.6 GeV theory error  *
*                                                                     *
*       other quark masses:                                           *
*                                                                     *
*   43) m_c(m_c)                                                      *
*   44) m_b(m_b)                                                      *
*                                                                     *
*       gauge couplings and related:                                  *
*                                                                     *
*   45) delta alpha_had^(3) (1.8 GeV)                                 *
*   46) (g_mu - 2 - alpha/pi)/2                                       *
*   47) tau lifetime                                                  *
*                                                                     *
*       neutrino-nucleon DIS:                                         *
*                                                                     *
*   48) g_L^2   (NuTeV 2002)                                          *
*   49) g_R^2   (NuTeV 2002)                                          *
*   50) kappa    (CCFR 1997)                                          *
*   51) R_nu    (CHARM 1984)                                          *
*   52) R_nu     (CDHS 1984)                                          *
*   53) R_nubar (CHARM 1984)                                          *
*   54) R_nubar  (CDHS 1984)                                          *
*   55) R_nubar  (CDHS 1979)                                          *
*                                                                     *
*       other low energy:                                             *
*                                                                     *
*   56) g_V^(nu,e) (CHARM II)                                         *
*   57) g_A^(nu,e) (CHARM II)                                         *
*   58) Q_W (e)    (SLAC E-158)                                       *
*   59) Q_W (e)    (e2ePV placeholder)                                *
*   60) Q_W (H)    (Qweak placeholder)                                *
*   61) Q_W (Cs)                                                      *
*   62) Q_W (Tl)                                                      *
*   63) 4*C1d + 9*C1u (polarized e- scattering; Young et al.)         *
*   64) 9*C1d - 4*C1d (polarized e- scattering; Young et al.)         *
*   65) CKM 1st row unitarity                                         * 
*   66) B(b --> s gamma)/B(b --> c e nu)                              *
*                                                                     *
*       other Tevatron:                                               *
*                                                                     *
*   67) A^FB (e) (CDF II)                                             *
*   68) PVDIS                                                         *
*                                                                     *
*       other LEP 2:                                                  *
*                                                                     *
*   69) sigma (hadrons) (183 GeV)                                     *
*   70) sigma (mu)      (183 GeV)                                     *
*   71) sigma (tau)     (183 GeV)                                     *
*   72) A^FB (mu)       (183 GeV)                                     *
*   73) A^FB (tau)      (183 GeV)                                     *
*   74) sigma (hadrons) (183 GeV)                                     *
*   75) sigma (mu)      (183 GeV)                                     *
*   76) sigma (tau)     (183 GeV)                                     *
*   77) A^FB (mu)       (183 GeV)                                     *
*   78) A^FB (tau)      (183 GeV)                                     *
*                                                                     *
*       new physics                                                   *
*                                                                     *
*   79) S parameter                                                   *
*   80) T parameter                                                   *
*                                                                     *
***********************************************************************

      data value/  91.1876d0,   2.4952d3,  41.5410d0,  20.8040d0,
     .             20.7850d0,  20.7640d0,   0.0145d0,   0.0169d0,
     .              0.0188d0,   0.1439d0,   0.1498d0,   0.0403d0,
     .              0.0980d0,   0.3710d0,  0.21629d0,   0.1721d0,
     .              0.0992d0,   0.0707d0,   0.9230d0,   0.6700d0,
     .             0.15138d0,   0.1544d0,   0.1420d0,   0.1360d0,
     .              0.1620d0,   0.8950d0,  80.3760d0,   2.1960d0,
     .             80.4320d0,   2.0570d0, 176.1000d0, 167.4000d0,
     .            186.0000d0, 180.1000d0, 168.4000d0, 170.9000d0,
     .            164.5000d0, 171.1000d0, 183.9000d0, 170.5000d0,
     .              0.1217d0, 172.4000d0,   0.0000d0,   0.0000d0,
     .             56.910d-4,  4511.07d0, 290.9300d0,   0.3010d0,
     .              0.0308d0,   0.5820d0,   0.3021d0,   0.3096d0,
c    .             0.03090d0,   0.5820d0,   0.3021d0,   0.3096d0,
     .              0.4030d0,   0.3840d0,   0.3650d0, - 0.0400d0,
c    .            - 0.5070d0, - 0.0403d0,   0.0713d0, -72.6200d0,
     .            - 0.5070d0, - 0.0403d0, - 0.0469d0,   0.0715d0,
     .            -73.1600d0,-116.4000d0, - 0.0285d0,   0.3420d0,
     .              1.0000d0, - 5.6400d0,  0.23158d0, - 0.8011d0,
     .             24.5400d0,   3.4400d0,   3.4300d0,   0.5470d0,
     .              0.6150d0,  22.3800d0,   3.1930d0,   3.1350d0,
     .              0.5620d0,   0.5970d0,   0.0000d0,   0.0000d0/

      data error/   0.0021d0,   0.0023d3,   0.0370d0,   0.0500d0,
     .              0.0330d0,   0.0450d0,   0.0025d0,   0.0013d0,
     .              0.0017d0,   0.0043d0,   0.0049d0,   0.0026d0,
     .              0.0110d0,   0.0220d0,   6.600d-4,   0.0030d0,
     .              0.0016d0,   0.0035d0,   0.0200d0,   0.0270d0,
     .              2.160d-3,   0.0060d0,   0.0150d0,   0.0150d0,
     .              0.0430d0,   0.0910d0,   0.0330d0,   0.0830d0,
     .              0.0390d0,   0.0620d0,   7.3000d0,  11.4000d0,
     .             11.5000d0,   5.3000d0,  12.8000d0,   2.5000d0,
     .              5.6000d0,   4.3000d0,  15.8000d0,   2.7000d0,
     .              1.000d-6,   1.3400d0,   0.0000d0,   0.0000d0,
     .              0.960d-4,   0.7400d0,   0.4800d0,   0.0015d0,  
     .              0.0011d0,   0.0041d0,   0.0041d0,   0.0043d0,
     .              0.0160d0,   0.0180d0,   0.0160d0,   0.0150d0,
c    .              0.0140d0,   0.1730d0,   0.0030d0,   0.4600d0,
     .              0.0140d0,   0.0053d0,   1.079d-3,   0.0029d0,
     .              0.3500d0,   3.6400d0,   0.0043d0,   0.0630d0,
     .              0.0006d0,   0.1400d0,   1.770d-3,   0.0048d0,
C Kai010            0.4300d0,   0.1400d0,   0.1800d0,   0.0000d0,  
     .              0.4300d0,   0.1400d0,   0.1800d0,   0.0340d0,  
     .              0.0440d0,   0.2500d0,   0.0830d0,   0.1020d0,
     .              0.0220d0,   0.0270d0,   1.000d-4,   0.0100d0/

      data ccls/ 
     .            0.01413d0, 0.04390d0, 0.06821d0,-0.09236d0,-0.00220d0,
     .           -0.00376d0,-0.02215d0,-0.04409d0,-0.03286d0,
     .            0.04390d0, 0.10309d0, 0.34215d0,-0.02890d0,-0.04929d0,
     .           -0.03337d0,-0.01819d0,-0.00435d0,-0.00296d0,
     .            0.06821d0, 0.34215d0, 0.14441d0,-0.12609d0,-0.13752d0,
     .           -0.09200d0,-0.04912d0,-0.00400d0,-0.00376d0,
     .           -0.09236d0,-0.02890d0,-0.12609d0, 0.18988d0,-0.06297d0,
     .           -0.03972d0, 0.44043d0,-0.00763d0,-0.00263d0,
     .           -0.00220d0,-0.04929d0,-0.13752d0,-0.06297d0, 0.02672d0,
     .           -0.05496d0,-0.02394d0,-0.01119d0, 0.00483d0,
     .           -0.00376d0,-0.03337d0,-0.09200d0,-0.03972d0,-0.05496d0,
     .            0.01442d0,-0.01761d0, 0.00094d0,-0.00882d0,
     .           -0.02215d0,-0.01819d0,-0.04912d0, 0.44043d0,-0.02394d0,
     .           -0.01761d0, 0.16417d0, 0.01983d0, 0.01762d0,
     .           -0.04409d0,-0.00435d0,-0.00400d0,-0.00763d0,-0.01119d0,
     .            0.00094d0, 0.01983d0, 0.00484d0,-0.04421d0,
     .           -0.03286d0,-0.00296d0,-0.00376d0,-0.00263d0, 0.00483d0,
     .           -0.00822d0, 0.01762d0,-0.04421d0, 0.00368d0/ 

       data cchf/
     . 0.05498d0, 0.17588d0, 0.10585d0,-0.07619d0, 0.07364d0,-0.03776d0,  
     . 0.17588d0, 0.04065d0,-0.03052d0, 0.05180d0,-0.03096d0, 0.05704d0,
     . 0.10585d0,-0.03052d0, 0.04011d0,-0.16622d0,-0.05564d0,-0.00370d0,
     .-0.07619d0, 0.05180d0,-0.16622d0, 0.03537d0, 0.02653d0,-0.03651d0,
     . 0.07364d0,-0.03096d0,-0.05564d0, 0.02653d0, 0.02397d0,-0.11795d0,
     .-0.03776d0, 0.05704d0,-0.00370d0,-0.03651d0,-0.11795d0, 0.01940d0/

       data ccslc/
     . 0.00721d0,-0.07381d0,-0.02898d0,-0.02920d0,
     .-0.07381d0, 0.00572d0,-0.01180d0,-0.00877d0,
     .-0.02898d0,-0.01180d0, 0.00108d0,-0.00601d0,
     .-0.02920d0,-0.00877d0,-0.00601d0, 0.00101d0/

       data ccdis/
     . 0.018d0,-0.059d0,-0.059d0,-0.059d0, 0.014d0, 0.014d0, 0.014d0,
     .-0.059d0, 0.307d0,-0.360d0,-0.360d0,-0.045d0,-0.045d0,-0.045d0,
     .-0.059d0,-0.360d0, 0.307d0,-0.360d0,-0.045d0,-0.045d0,-0.045d0,
     .-0.059d0,-0.360d0,-0.360d0, 0.307d0,-0.045d0,-0.045d0,-0.045d0,
     . 0.014d0,-0.045d0,-0.045d0,-0.045d0, 0.051d0,-0.125d0,-0.125d0,
     . 0.014d0,-0.045d0,-0.045d0,-0.045d0,-0.125d0, 0.051d0,-0.125d0,
     . 0.014d0,-0.045d0,-0.045d0,-0.045d0,-0.125d0,-0.125d0, 0.051d0/

       data ccmt/
     . 0.351d0,-0.212d0,-0.262d0,-0.159d0,-0.031d0,-0.122d0,-0.270d0,
     .-0.041d0,-0.056d0,-0.049d0,-0.104d0,
     .-0.212d0, 0.136d0,-0.090d0,-0.043d0,-0.027d0,-0.062d0,-0.118d0,
     .-0.048d0,-0.006d0,-0.009d0,-0.051d0,
     .-0.262d0,-0.090d0, 0.161d0,-0.024d0,-0.014d0, 0.006d0,-0.156d0,
     .-0.057d0, 0.006d0, 0.012d0,-0.061d0,
     .-0.159d0,-0.043d0,-0.024d0, 0.159d0,-0.122d0,-0.076d0,-0.120d0,
     .-0.024d0,-0.033d0,-0.020d0,-0.174d0,
     .-0.031d0,-0.027d0,-0.014d0,-0.122d0, 0.041d0,-0.030d0,-0.053d0,
     .-0.022d0, 0.003d0,-0.001d0,-0.049d0,
     .-0.122d0,-0.062d0, 0.006d0,-0.076d0,-0.030d0, 0.095d0,-0.048d0,
     .-0.141d0,-0.028d0,-0.137d0, 0.021d0,
     .-0.270d0,-0.118d0,-0.156d0,-0.120d0,-0.053d0,-0.048d0, 0.295d0,
     .-0.042d0,-0.003d0, 0.068d0,-0.277d0,
     .-0.041d0,-0.048d0,-0.057d0,-0.024d0,-0.022d0,-0.141d0,-0.042d0,
     . 0.052d0,-0.025d0,-0.061d0, 0.026d0,
     .-0.056d0,-0.006d0, 0.006d0,-0.033d0, 0.003d0,-0.028d0,-0.003d0,
     .-0.025d0, 0.009d0,-0.030d0, 0.024d0,
     .-0.049d0,-0.009d0, 0.012d0,-0.020d0,-0.001d0,-0.137d0, 0.068d0,
     .-0.061d0,-0.030d0, 0.104d0,-0.277d0,
     .-0.104d0,-0.051d0,-0.061d0,-0.174d0,-0.049d0, 0.021d0,-0.277d0,
     . 0.026d0, 0.024d0,-0.277d0, 0.250d0/

       cctau =   0.012d0
       ccmgw = - 0.177d0
       cchad = - 0.170d0 ! - sign: prediction error is added to exp. error
       ccntv = - 0.017d0
       ccnue = - 0.050d0

*   fit parameters:
                                 mz = xval(1)
      if (flagmt.eqv..true.)     mt = xval(2)
                                 mb = xval(3)
      if (flagmc.eqv..true.)     mc = xval(4)
                             alfas0 = xval(5)
      if (flagal.eqv..true.) dahad3 = xval(6)
      if (flagmh.eqv..true.)     mh = dexp(xval(7))
      if (flagT.eqv..true.)    Tpar = xval(8)
      if (flagS.eqv..true.)    Spar = xval(9) 
                               Upar = xval(10)
      if (flgrho.eqv..true.)   Brho = xval(11)
      if (fkappa.eqv..true.) Bkappa = xval(12)
                               Zpar = xval(15)
C     if (fzprim.eqv..true.)    mzp = dexp(xval(16))
      if (fzprim.eqv..true.)    mzp = 1.d3/xval(16)
      if (fsinth.eqv..true.)  sinth = xval(17)
                             lambdg = xval(18)

C Kai008 -------------------------------------------------
C     Initializes the values of the input parameters.
C --------------------------------------------------------
      if (flgfitx.eqv..true.)  fitx = dexp(xval(28))
      if (flgtph.eqv..true.) fittph = xval(29)
      if (flgs2b.eqv..true.) fits2b = xval(30)

        fitsph = fittph / (1.0d0+fittph)
        fitcph = 1.0d0 - fitsph

        kkem   = 1.d0/127.918
        kkcc   = 1.d0/2.d0*
     .  (1.d0 +dsqrt(1.d0-4.d0*pi1*kkem/dsqrt(2.d0)/mz**2/gf))
        kkss   = 1.d0 - kkcc
C ---------------------------------------------------------

c       print*,"xval=",mz,mt,mb,mc,alfas0,dahad3,mh,Tpar,Spar,Upar,Brho,
c     $  Bkappa, Zpar, mzp, sinth,lambdg
c       print*,"---------------------------------------"
c        print*,"logical=",flagmt,flagmc,flagal,flgfitx,flgtph,
c     $  flgs2b,flagmh
c        print*,"---------------------------------------"
C   The scale mu0 at which alphas is initialized (alfas0) 
C   MUST be below mt and MUST NOT be below mb.

       SWpar = Spar + Upar
            
       mz2 = mz**2
      mzp2 = mzp**2
       mt2 = mt**2
       mh2 = mh**2
       mu0 = mz

      rathz2 = mh2/mz2
      rattz2 = mt2/mz2
      ratth2 = mt2/mh2      

C   For specified Higgs sector in Z prime models (flgzph = .true.):

      flgzph = .false.
      ratgRL = lambdg
      sigma  = 1.d0

      if (flagzp.ne.0) then
         if (flgzph.eqv..true.) then
            if (flagzp.eq.1) sinth = (1 - 3*sigma/2)/(1 + sigma)
     .           *dsqrt(0.4626d0*lambdg/3)*mz2/mzp2
            if (flagzp.eq.2) sinth = (1 - sigma)/(1 + sigma)
     .           *dsqrt(0.257d0*lambdg)*mz2/mzp2
            if (flagzp.eq.3) sinth = (4*sigma - 1)/(sigma + 1)
     .           *dsqrt(0.0257d0*lambdg)*mz2/mzp2
            if (flagzp.eq.4) sinth = (1 - sigma/(ratgRL**2*3.3234d0-1))/
     .           (1+sigma)*dsqrt(ratgRL**2*0.7687d0 - 0.2313d0)*mz2/mzp2
            if (flagzp.eq.13) sinth = 0.d0
         endif
         sinth2 = sinth**2
         costh2 = 1.d0 - sinth2
          costh = dsqrt(costh2)
           mz02 = costh2*mz2 + sinth2*mzp2            
      endif

      call sin2thetaw
 
      if (flagzp.ne.0) then
         if (flagzp.eq.4) then
            ratg21 = dsqrt(5*sinhat/3)
         else if (flagzp.eq.6) then
            ratg21 = 0.20d0
         else if (flagzp.eq.13) then
C           ratg21 = dsqrt(sinhat*coshat/4/pi1/alphat)
            ratg21 = dsqrt(0.23112d0*0.76888d0/4/pi1*127.92d0)
         else
            ratg21 = dsqrt(5*sinhat*lambdg/3)
         endif
         if (flgzph.eqv..true.) then
            if (flagzp.eq.1) chiggs = 
     .           (1 - 3*sigma/2)/(1 + sigma)*2/dsqrt(10.d0)
            if (flagzp.eq.2) chiggs = 
     .           (1 - sigma)/(1 + sigma)*dsqrt(2/3.d0)
            if (flagzp.eq.3) chiggs = 
     .           (4*sigma - 1)/(sigma + 1)/dsqrt(15.d0)
            do 150 i = 1,2
               if (flagzp.eq.4) chiggs = 1.d0/(1 + sigma)
     .              *(1 - sigma/(ratgRL**2*coshat/sinhat - 1))
     .              *dsqrt(3*(ratgRL**2*coshat/sinhat - 1)/5)
               
               sinth = chiggs*ratg21*mz2/mzp2
              sinth2 = sinth**2
              costh2 = 1.d0 - sinth2
               costh = dsqrt(costh2)
                mz02 = costh2*mz2 + sinth2*mzp2  
  
                call sin2thetaw

                if (flagzp.eq.4) then
                   ratg21 = dsqrt(5*sinhat/3)
                else if (flagzp.eq.6) then
                   ratg21 = 0.20d0
                else if (flagzp.eq.13) then
                   ratg21 = dsqrt(sinhat*coshat/4/pi1/alpha)
                else
                   ratg21 = dsqrt(5*sinhat*lambdg/3)
                endif
 150        continue             
         endif
         rhoeff = mz02*          (costh2/mz2 + sinth2/mzp2)
         rhoezp = mz02*ratg21**2*(sinth2/mz2 + costh2/mzp2)
         rhozzp = mz02*ratg21*sinth*costh*(1/mz2 - 1/mzp2)
         if (fobliq.eqv..true.) then
            rhoeff = rhoeff/(1 - alphat*Tpar)
            rhoezp = rhoezp/(1 - alphat*Tpar)
            rhozzp = rhozzp/(1 - alphat*Tpar)
         endif
C        xval(19) = dsin(ima*pi1/ssteps)  ! for plot of model 13
C        xval(20) = dcos(ima*pi1/ssteps)  ! for plot of model 13
         call zprimecoup(xval,npar)
      endif

      call z0pole(gamma,sigmah,R,sin2te,alr,afb,0,9)
      call rho0
      call sin2theta0
      call rhof(rho,0,9)
      call kappaf(kappa,0,9)
      call wwprod(gammaw)
      call lep200(xs200,afb200)
      call smrule(sumrlc,momntc,errorc,sumrlb,momntb,errorb,ccmbc)

      value(43) = sumrlc(2)
      value(44) = sumrlb(6)
      error(43) = errorc(2)
      error(44) = errorb(6)

      do 200 i = 1, 80
         smval(i) = 0.d0
          pull(i) = 0.d0
 200  continue

      smval( 1) =     mz
      smval( 2) =  gamma(11)
      smval( 3) = sigmah
      smval( 4) =      R(1)
      smval( 5) =      R(2)
      smval( 6) =      R(3)
      smval( 7) =    afb(1)
      smval( 8) =    afb(2)
      smval( 9) =    afb(3)

      smval(10) =    alr(3)
      smval(11) =    alr(1)
      smval(12) =    afb(10)
      smval(13) =    afb(8)
      smval(14) =      R(8)/(R(4) + R(7) + R(8))

      smval(15) =      R(9)
      smval(16) =      R(5)
      smval(17) =    afb(9)
      smval(18) =    afb(5)
      smval(19) =    alr(9)
      smval(20) =    alr(5)

      smval(21) =    alr(1)
      smval(22) =    alr(1)
      smval(23) =    alr(2)
      smval(24) =    alr(3)
      smval(25) =    alr(1)
      smval(26) =    alr(8) 

      smval(27) =     mw
      smval(28) = gammaw(7)
      smval(29) =     mw
      smval(30) = gammaw(7)

      call polemasses(6,smval(31))
      call polemasses(6,smval(32))
      call polemasses(6,smval(33))
      call polemasses(6,smval(34))
      call polemasses(6,smval(35))
      call polemasses(6,smval(36))
      call polemasses(6,smval(37))
      call polemasses(6,smval(38))
      call polemasses(6,smval(39))
      call polemasses(6,smval(40))
C     call polemasses(6,smval(41))
      smval(41) = alfas0 
      call polemasses(6,smval(42))

      smval(43) = momntc(2)
      smval(44) = momntb(6)

      smval(45) = dahad3
      call anomagmntmu(smval(46))
      call taulifetime(smval(47),therr2)

      call nuhnutev(dummy,epsu_L,epsd_L,epsu_R,epsd_R)
      smval(48) = epsu_L**2 + epsd_L**2
      smval(49) = epsu_R**2 + epsd_R**2
      call nuhccfr (smval(50))
      call nuhcdhs (smval(51),smval(53),dummy)
      call nuhcdhs (smval(52),smval(54),smval(55))

      call nue(0.d0,smval(56),smval(57))
      call moller(smval(58),0.0260d0,0.599d0)
      call moller(smval(59),0.0056d0,0.571d0)
      call apv(smval(60),1,1,C1u,C1d,C2u,C2d)
      call apv(smval(61),55,133,C1u,C1d,C2u,C2d)
      call apv(smval(62),81,205,C1u,C1d,C2u,C2d)
      smval(63) =   0.9137*C1u + 0.4065*C1d ! email Ross Young 11/13/08
      smval(64) = - 0.4065*C1u + 0.9137*C1d ! email Ross Young 11/13/08

      smval(65) = 1.d0
      if (flagzp.ne.0) smval(65) = 1.d0 + 
     .     3*ratg21**2*alphat/pi1/sinhat/coshat*eps2_L(1)*
     .     (eps2_L(1) - eps2_L(4))*dlog(mzp2/mw2)/(mzp2/mw2 - 1)

      call bsgamma(sgcenu)
      smval(66) = dlog(sgcenu)

      smval(67) = sin2te(1)
      smval(68) = (2*C1u - C1d) + 0.84d0*(2*C2u - C2d)
      smval(68) = (2*C1u - C1d) + 0.84d0*(2*C2u - C2d)

C     smval(69) =  xs200(1,10)
C     smval(70) =  xs200(1, 2)
C     smval(71) =  xs200(1, 3)
C     smval(72) = afb200(1, 2)
C     smval(73) = afb200(1, 3)
C     smval(74) =  xs200(2,10)
C     smval(75) =  xs200(2, 2)
C     smval(76) =  xs200(2, 3)
C     smval(77) = afb200(2, 2)
C     smval(78) = afb200(2, 3)

      smval(79) = Spar
      smval(80) = Tpar

      chi2 = 0.d0
c      print*,"___________________________"
c      print*,"chi2=",chi2,sigmah,R
c      write(*,*)("smval[",ii2,"]=",smval(ii2),ii2=1,68)
C Kai009 ------------------------------
C    Changed set of observables that are used
C    for the global fit.
C -------------------------------------

C    Z lineshape, Other LEP 1 measurements, LEP 1 and SLC heavy flavor,
C    Other SLD asymmetries.
      do 300 i = 2,26
         prtpll(i) = 1
         pull(i) = (value(i) - smval(i))/error(i)
         chi2 = chi2 + pull(i)**2
 300  continue

C    W properties
      do 302 i = 27,30
         prtpll(i) = 1
         pull(i) = (value(i) - smval(i))/error(i)
         chi2 = chi2 + pull(i)**2
 302  continue
     
C    m_t (total TeVatron average), Other quark masses, Gauge couplings
C    and related (tau lifetime excluded)
C     do 303 i = 42,46
      do 303 i = 42,42
         prtpll(i) = 1
         pull(i) = (value(i) - smval(i))/error(i)
         chi2 = chi2 + pull(i)**2
 303  continue
   
C    tau lifetime 
      prtpll(47) = 1
      pull(47) = (value(47) - smval(47))/dsqrt((error(47)**2 + therr2))
      chi2 = chi2 + pull(47)**2

C    Neutrino-nucleon DIS, Other low energy (first three)
C     do 304 i = 48,58
C        prtpll(i) = 1
C        pull(i) = (value(i) - smval(i))/error(i)
C        chi2 = chi2 + pull(i)**2
C304  continue

      do 304 i = 48,58
         prtpll(i) = 1
         pull(i) = (value(i) - smval(i))/error(i)
         chi2 = chi2 + pull(i)**2
 304  continue

C    Upcoming measurements of weak charges (e and p)
      do 305 i = 59,60
         prtpll(i) = 1
         pull(i) = (value(i) - smval(i))/error(i)
         chi2 = chi2 + pull(i)**2
 305  continue

C    Other low energy (last six), A^FB(e) from TeVatron (CDF II)
      do 306 i = 61,64
          prtpll(i) = 1
          pull(i) = (value(i) - smval(i))/error(i)
          chi2 = chi2 + pull(i)**2
 306  continue

C    Other LEP 2
C     do 307 i = 69,78
C         prtpll(i) = 1
C         pull(i) = (value(i) - smval(i))/error(i)
C         chi2 = chi2 + pull(i)**2
C307  continue

C Kai009 ---------------------------------------

C Kai022 ---------------------------------------
C     Write initial pulls in the *.pull file
C ----------------------------------------------
      if (flfout.eqv..true.) then   
         if (plltr.eq.1) then
        
            do 800 i = 1,80
               if (prtpll(i).eq.1) then
                  write(50,FMT='(I3,A,F8.4)') i,'\t', pull(i)
               endif
 800        continue
      
            plltr = 0

         endif
      endif

C Kai022 ---------------------------------------


       flcorr = .true.

      if (flcorr.eqv..true.) then
         do 410 i = 1, 9
            chi2 = chi2 + ccls(i,i)*pull(i)**2
            do 400 j = i+1, 9
               chi2 = chi2 + 2*ccls(j,i)*pull(i)*pull(j)
 400        continue
 410     continue
         
         chi2 = chi2 - pull(10)**2 - pull(11)**2 + 1/dsqrt(1 - cctau**2)
     .        *(pull(10)**2 + pull(11)**2 - 2*cctau*pull(10)*pull(11))

         do 430 i = 15, 20
            chi2 = chi2 + cchf(i,i)*pull(i)**2
            do 420 j = i+1, 20
               chi2 = chi2 + 2*cchf(i,j)*pull(i)*pull(j)
 420        continue
 430     continue
         
         do 450 i = 21, 24
            chi2 = chi2 + ccslc(i,i)*pull(i)**2
            do 440 j = i+1, 24
               chi2 = chi2 + 2*ccslc(i,j)*pull(i)*pull(j)
 440        continue
 450     continue

ccccc         chi2 = chi2 - pull(29)**2 - pull(30)**2 + 1/dsqrt(1 - ccmgw**2)
ccccc     .        *(pull(29)**2 + pull(30)**2 - 2*ccmgw*pull(29)*pull(30))

ccccc        chi2 = chi2 + ccmt(i,i)*pull(i)**2 
ccccc        do 460 j = i+1, 41
ccccc           chi2 = chi2 + 2*ccmt(j,i)*pull(i)*pull(j)
c 460        continue
c 470     continue

         chi2 = chi2 - pull(43)**2 - pull(44)**2 + 1/dsqrt(1 - ccmbc**2)
     .        *(pull(43)**2 + pull(44)**2 - 2*ccmbc*pull(43)*pull(44))

         chi2 = chi2 - pull(45)**2 - pull(46)**2 + 1/dsqrt(1 - cchad**2)
     .        *(pull(45)**2 + pull(46)**2 - 2*cchad*pull(45)*pull(46))

         chi2 = chi2 - pull(48)**2 - pull(49)**2 + 1/dsqrt(1 - ccntv**2)
     .        *(pull(48)**2 + pull(49)**2 - 2*ccntv*pull(48)*pull(49))

         do 490 i = 50, 55
            chi2 = chi2 + ccdis(i,i)*pull(i)**2
            do 480 j = i+1, 55
               chi2 = chi2 + 2*ccdis(i,j)*pull(i)*pull(j)
 480        continue
 490     continue

         chi2 = chi2 - pull(56)**2 - pull(57)**2 + 1/dsqrt(1 - ccnue**2)
     .        *(pull(56)**2 + pull(57)**2 - 2*ccnue*pull(56)*pull(57))
         
      endif
      
      if (fhiggs.eqv..true.) then

C  Tevatron Higgs searches: FERMILAB-PUB-08-068-E (M_H < 155 GeV) and
C                           FERMILAB-PUB-08-270-E (M_H > 155 GeV)

         if (mh.lt.205.d0.and.mh.ge.200.d0)
     .                  chi2 = chi2 + (205.d0 - mh)*0.02d0
         if (mh.lt.200.d0.and.mh.ge.195.d0)
     .                  chi2 = chi2                              + 0.1d0
         if (mh.lt.195.d0.and.mh.ge.190.d0) 
     .                  chi2 = chi2 + (195.d0 - mh)*0.14d0       + 0.1d0
         if (mh.lt.190.d0.and.mh.ge.185.d0)
     .                  chi2 = chi2 + (190.d0 - mh)*0.16d0       + 0.8d0
         if (mh.lt.185.d0.and.mh.ge.180.d0) 
     .                  chi2 = chi2 + (185.d0 - mh)*0.38d0       + 1.6d0
         if (mh.lt.180.d0.and.mh.ge.175.d0) 
     .                  chi2 = chi2 + (180.d0 - mh)*0.18d0       + 3.5d0
         if (mh.lt.175.d0.and.mh.ge.170.d0) 
     .                  chi2 = chi2 + (175.d0 - mh)*0.28d0       + 4.4d0
         if (mh.lt.170.d0.and.mh.ge.165.d0) 
     .                  chi2 = chi2 - (170.d0 - mh)*0.38d0       + 5.8d0
         if (mh.lt.165.d0.and.mh.ge.160.d0) 
     .                  chi2 = chi2 - (165.d0 - mh)*0.28d0       + 3.9d0
         if (mh.lt.160.d0.and.mh.ge.155.d0) 
     .                  chi2 = chi2 - (160.d0 - mh)*0.08d0       + 2.5d0
         if (mh.lt.155.d0.and.mh.ge.150.d0) 
     .                  chi2 = chi2 + (155.d0 - mh)*0.74d0/3.d0  + 2.1d0
         if (mh.lt.150.d0.and.mh.ge.140.d0) 
     .                  chi2 = chi2 - (150.d0 - mh)*70/240.d0 + 80/24.d0
         if (mh.lt.140.d0.and.mh.ge.130.d0) 
     .                  chi2 = chi2 - (140.d0 - mh)*20/240.d0 + 10/24.d0
         if (mh.lt.130.d0.and.mh.ge.120.d0) 
     .                  chi2 = chi2 - (130.d0 - mh)*8/240.d0  - 10/24.d0
         if (mh.lt.120.d0.and.mh.ge.115.d0) 
     .                  chi2 = chi2 + (120.d0 - mh)*52/240.d0 - 18/24.d0
         if (mh.lt.115.d0) 
     .                  chi2 = chi2 + (115.d0 - mh)*22/240.d0 +  8/24.d0

C  final LEP 2 Higgs searches:

         if (mh.lt.120.0d0.and.mh.ge.119.6d0) 
     .                   chi2 = chi2 + (120.0d0 - mh)*0.5d0
         if (mh.lt.119.6d0.and.mh.ge.118.8d0) 
     .                   chi2 = chi2 - (119.6d0 - mh)         +  0.2d0
         if (mh.lt.118.8d0.and.mh.ge.118.2d0) 
     .                   chi2 = chi2 - (118.8d0 - mh)/3.d0    -  0.6d0
         if (mh.lt.118.2d0.and.mh.ge.117.0d0) 
     .                   chi2 = chi2 - (118.2d0 - mh)*0.75d0  -  0.8d0
         if (mh.lt.117.0d0.and.mh.ge.116.4d0) 
     .                   chi2 = chi2 + (117.0d0 - mh)*2/3.d0  -  1.7d0
         if (mh.lt.116.4d0.and.mh.ge.115.8d0) 
     .                   chi2 = chi2 - (116.4d0 - mh)/3.d0    -  1.3d0
         if (mh.lt.115.8d0.and.mh.ge.115.2d0) 
     .                   chi2 = chi2 + (115.8d0 - mh)*2/3.d0  -  1.5d0
         if (mh.lt.115.2d0.and.mh.ge.114.6d0) 
     .                   chi2 = chi2 + (115.2d0 - mh)*10/3.d0 -  1.1d0
         if (mh.lt.114.6d0.and.mh.ge.113.4d0) 
     .                   chi2 = chi2 + (114.6d0 - mh)*6.d0    +  0.9d0
         if (mh.lt.113.4d0.and.mh.ge.112.2d0) 
     .                   chi2 = chi2 + (113.4d0 - mh)*23/3.d0 +  8.1d0
         if (mh.lt.112.2d0.and.mh.ge.111.6d0) 
     .                   chi2 = chi2 + (112.2d0 - mh)*16/3.d0 + 17.3d0
         if (mh.lt.111.6d0.and.mh.ge.111.0d0) 
     .                   chi2 = chi2 + (111.6d0 - mh)*41/6.d0 + 20.5d0
         if (mh.lt.111.0d0.and.mh.ge.110.4d0) 
     .                   chi2 = chi2 + (111.0d0 - mh)*9.5d0   + 24.6d0
         if (mh.lt.110.4d0.and.mh.ge.109.8d0) 
     .                   chi2 = chi2 + (110.4d0 - mh)*31/3.d0 + 30.3d0
         if (mh.lt.109.8d0.and.mh.ge.109.2d0)
     .                   chi2 = chi2 + (109.8d0 - mh)*6.5d0   + 36.5d0
         if (mh.lt.109.2d0.and.mh.ge.108.6d0) 
     .                   chi2 = chi2 + (109.2d0 - mh)*6.d0    + 40.4d0
         if (mh.lt.108.6d0.and.mh.ge.108.0d0) 
     .                   chi2 = chi2 + (108.6d0 - mh)*43/6.d0 + 44.0d0
         if (mh.lt.108.0d0) chi2 = chi2 + (108.0d0 - mh)*8.5d0+ 48.3d0

      endif

      return
      end

***************************************************************************
*                                                                         *
*                 Earlier exclusion and search results                    *
*                                                                         *
***************************************************************************

C  Tevatron Higgs searches: FERMILAB-PUB-08-068-E
      
C         if (mh.lt.201.d0.and.mh.ge.190.d0) 
C     .                  chi2 = chi2 + (200.d0 - mh)*30/240.d0 +  3/24.d0
C         if (mh.lt.190.d0.and.mh.ge.180.d0) 
C     .                  chi2 = chi2 - (190.d0 - mh)/240.d0    + 33/24.d0
C         if (mh.lt.180.d0.and.mh.ge.170.d0) 
C     .                  chi2 = chi2 + (180.d0 - mh)*51/240.d0 + 32/24.d0
C         if (mh.lt.170.d0.and.mh.ge.160.d0) 
C     .                  chi2 = chi2 + (170.d0 - mh)*14/240.d0 + 83/24.d0
C         if (mh.lt.160.d0.and.mh.ge.150.d0) 
C     .                  chi2 = chi2 - (160.d0 - mh)*17/240.d0 + 97/24.d0
C         if (mh.lt.150.d0.and.mh.ge.140.d0) 
C     .                  chi2 = chi2 - (150.d0 - mh)*70/240.d0 + 80/24.d0
C         if (mh.lt.140.d0.and.mh.ge.130.d0) 
C     .                  chi2 = chi2 - (140.d0 - mh)*20/240.d0 + 10/24.d0
C         if (mh.lt.130.d0.and.mh.ge.120.d0) 
C     .                  chi2 = chi2 - (130.d0 - mh)*8/240.d0  - 10/24.d0
C         if (mh.lt.120.d0.and.mh.ge.115.d0) 
C     .                  chi2 = chi2 + (120.d0 - mh)*52/240.d0 - 18/24.d0
C         if (mh.lt.115.d0) 
C     .                  chi2 = chi2 + (115.d0 - mh)*22/240.d0 +  8/24.d0

C  summer 2001 Higgs searches (LEP):

C      if (mh.lt.120.5d0.and.mh.ge.117.5d0) 
C     .                   chi2 = chi2 -  (120.5d0 - mh)*0.7d0
C      if (mh.lt.117.5d0.and.mh.ge.116.5d0) 
C     .                   chi2 = chi2                         - 2.1d0
C      if (mh.lt.116.5d0.and.mh.ge.116.0d0) 
C     .                   chi2 = chi2 -  (116.5d0 - mh)*1.8d0 - 2.1d0
C      if (mh.lt.116.0d0.and.mh.ge.115.0d0) 
C     .                   chi2 = chi2                         - 3.0d0
C      if (mh.lt.115.0d0.and.mh.ge.114.0d0) 
C     .                   chi2 = chi2 +  (115.0d0 - mh)*3.8d0 - 3.0d0
C      if (mh.lt.114.0d0.and.mh.ge.112.0d0) 
C     .                   chi2 = chi2 +  (114.0d0 - mh)*6.0d0 + 0.8d0
C      if (mh.lt.112.0d0) chi2 = chi2 +  (112.0d0 - mh)*8.1d0 +12.8d0

C  winter 2001 Higgs searches (LEP):

C      if (mh.lt.121.5d0.and.mh.ge.116.5d0) 
C     .                   chi2 = chi2 -  (121.5d0 - mh)*1.0d0
C      if (mh.lt.116.5d0.and.mh.ge.115.0d0) 
C     .                   chi2 = chi2 -  (116.5d0 - mh)*1.2d0 - 5.0d0
C      if (mh.lt.115.0d0.and.mh.ge.113.5d0) 
C     .                   chi2 = chi2 +  (115.0d0 - mh)*3.0d0 - 6.8d0
C      if (mh.lt.113.5d0.and.mh.ge.110.5d0) 
C     .                   chi2 = chi2 +  (113.5d0 - mh)*5.8d0 - 2.3d0
C      if (mh.lt.110.5d0.and.mh.ge.110.0d0) 
C     .                   chi2 = chi2 +  (110.5d0 - mh)*13.4d0+15.1d0
C      if (mh.lt.110.0d0.and.mh.ge.109.2d0) 
C     .                   chi2 = chi2 +  (110.0d0 - mh)*0.5d0 +21.8d0
C      if (mh.lt.109.2d0) chi2 = chi2 +  (109.2d0 - mh)*8.1d0 +22.2d0

C  fall 2000 Higgs searches (LEP):

C      if (mh.lt.122.0d0.and.mh.ge.118.5d0) 
C     .                   chi2 = chi2 -  (122.0d0 - mh)*0.6d0
C      if (mh.lt.118.5d0.and.mh.ge.117.5d0) 
C     .                   chi2 = chi2 -  (118.5d0 - mh)*0.5d0 - 2.1d0
C      if (mh.lt.117.5d0.and.mh.ge.116.5d0) 
C     .                   chi2 = chi2 -  (117.5d0 - mh)*0.6d0 - 2.6d0
C      if (mh.lt.116.5d0.and.mh.ge.115.5d0) 
C     .                   chi2 = chi2 -  (116.5d0 - mh)*1.4d0 - 3.2d0
C      if (mh.lt.115.5d0.and.mh.ge.115.0d0) 
C     .                   chi2 = chi2 -  (115.5d0 - mh)*3.2d0 - 4.6d0
C      if (mh.lt.115.0d0.and.mh.ge.114.5d0) 
C     .                   chi2 = chi2 +  (115.0d0 - mh)*3.4d0 - 6.2d0
C      if (mh.lt.114.5d0.and.mh.ge.114.0d0) 
C     .                   chi2 = chi2 +  (114.5d0 - mh)*1.6d0 - 4.5d0
C      if (mh.lt.114.0d0.and.mh.ge.112.5d0) 
C     .                   chi2 = chi2 +  (114.0d0 - mh)*4.2d0 - 3.7d0
C      if (mh.lt.112.5d0.and.mh.ge.111.5d0) 
C     .                   chi2 = chi2 +  (112.5d0 - mh)*6.3d0 + 2.6d0
C      if (mh.lt.111.5d0.and.mh.ge.111.0d0) 
C     .                   chi2 = chi2 +  (111.5d0 - mh)*1.0d0 + 8.9d0
C      if (mh.lt.111.0d0.and.mh.ge.110.5d0) 
C     .                   chi2 = chi2 +  (111.0d0 - mh)*9.4d0 + 9.4d0
C      if (mh.lt.110.5d0.and.mh.ge.110.0d0) 
C     .                   chi2 = chi2 +  (110.5d0 - mh)*2.8d0 +14.1d0
C      if (mh.lt.110.0d0.and.mh.ge.109.5d0) 
C     .                   chi2 = chi2 +  (110.0d0 - mh)*7.4d0 +15.5d0
C      if (mh.lt.109.5d0.and.mh.ge.109.0d0) 
C     .                   chi2 = chi2 +  (109.5d0 - mh)*2.2d0 +19.2d0
C      if (mh.lt.109.0d0) chi2 = chi2 +  (109.0d0 - mh)*7.8d0 +20.3d0

C  summer 2000 Higgs searches (LEP):

C      if (mh.lt.123.2d0.and.mh.ge.118.2d0) 
C     .                   chi2 = chi2 -  (123.2d0 - mh)*0.4d0
C      if (mh.lt.118.2d0.and.mh.ge.115.2d0) 
C     .                   chi2 = chi2 -  (118.2d0 - mh)*1.2d0 - 2.0d0
C      if (mh.lt.115.2d0.and.mh.ge.114.4d0) 
C     .                   chi2 = chi2                         - 5.6d0
C      if (mh.lt.114.4d0.and.mh.ge.113.4d0) 
C     .                   chi2 = chi2 +  (114.4d0 - mh)*1.6d0 - 5.6d0
C      if (mh.lt.113.4d0.and.mh.ge.111.4d0) 
C     .                   chi2 = chi2 +  (113.4d0 - mh)*3.3d0 - 4.0d0
C      if (mh.lt.111.4d0.and.mh.ge.110.4d0) 
C     .                   chi2 = chi2 +  (111.4d0 - mh)*1.0d0 + 2.6d0
C      if (mh.lt.110.4d0) chi2 = chi2 +  (110.4d0 - mh)*5.6d0 + 3.6d0

C  fall 1999 Higgs exclusion (LEP):

C      fac = 2*dlog(10.d0)/100
C      if (mh.lt.111.5d0.and.mh.ge.106.0d0) 
C     .                   chi2 = chi2 +  (111.5d0 - mh)*4*fac
C      if (mh.lt.106.0d0.and.mh.ge.104.0d0) 
C     .                   chi2 = chi2 + ((106.0d0 - mh)*19 +  22)*fac
C      if (mh.lt.104.0d0.and.mh.ge.103.0d0) 
C     .                   chi2 = chi2 + ((104.0d0 - mh)*70 +  60)*fac
C      if (mh.lt.103.0d0.and.mh.ge.101.0d0) 
C     .                   chi2 = chi2 + ((103.0d0 - mh)*45 + 130)*fac
C      if (mh.lt.101.0d0.and.mh.ge.99.0d0) 
C     .                   chi2 = chi2 + ((101.0d0 - mh)*10 + 220)*fac
C      if (mh.lt.99.0d0.and.mh.ge.98.5d0) 
C     .                   chi2 = chi2 +  ((99.0d0 - mh)*140+ 240)*fac
C      if (mh.lt.98.5d0.and.mh.ge.97.5d0) 
C     .                   chi2 = chi2 +                      310 *fac
C      if (mh.lt.97.5d0.and.mh.ge.97.0d0) 
C     .                   chi2 = chi2 +  ((97.5d0 - mh)*220+ 310)*fac
C      if (mh.lt.97.0d0.and.mh.ge.96.5d0) 
C     .                   chi2 = chi2 +  ((97.0d0 - mh)*60 + 420)*fac
C      if (mh.lt.96.5d0.and.mh.ge.95.5d0) 
C     .                   chi2 = chi2 +  ((96.5d0 - mh)*150+ 450)*fac
C      if (mh.lt.95.5d0.and.mh.ge.94.5d0) 
C     .                   chi2 = chi2 +  ((95.5d0 - mh)*70 + 600)*fac
C      if (mh.lt.94.5d0)  chi2 = chi2 +  ((94.5d0 - mh)*180+ 670)*fac

C  summer 1999 Higgs exclusion (LEP):

C      fac = 2*dlog(10.d0)/100
C      if (mh.lt.98.0d0.and.mh.ge.97.0d0) 
C     .                   chi2 = chi2 + ( 98.0d0 - mh)*22*fac
C      if (mh.lt.97.0d0.and.mh.ge.96.0d0) 
C     .                   chi2 = chi2 + ((97.0d0 - mh)*50 +  22)*fac
C      if (mh.lt.96.0d0.and.mh.ge.94.0d0) 
C     .                   chi2 = chi2 + ((96.0d0 - mh)*79 +  72)*fac
C      if (mh.lt.94.0d0.and.mh.ge.93.0d0) 
C     .                   chi2 = chi2 + ((94.0d0 - mh)*130+ 230)*fac
C      if (mh.lt.93.0d0.and.mh.ge.92.0d0) 
C     .                   chi2 = chi2 + ((93.0d0 - mh)*120+ 360)*fac
C      if (mh.lt.92.0d0)  chi2 = chi2 + ((92.0d0 - mh)*196+ 480)*fac

C  winter 1999 Higgs exclusion (L3):

C      fac = 2*dlog(10.d0)/100
C      if (mh.lt.101.5d0.and.mh.ge.97.5d0) 
C     .                   chi2 = chi2 + (101.5d0 - mh)*13*fac
C      if (mh.lt.97.5d0.and.mh.ge.96.0d0) 
C     .                   chi2 = chi2 + ((97.5d0 - mh)*32 +  52)*fac
C      if (mh.lt.96.0d0.and.mh.ge.94.0d0) 
C     .                   chi2 = chi2 + ((96.0d0 - mh)*40 + 100)*fac
C      if (mh.lt.94.0d0.and.mh.ge.93.0d0) 
C     .                   chi2 = chi2 + ((94.0d0 - mh)*20 + 180)*fac
C      if (mh.lt.93.0d0.and.mh.ge.92.0d0) 
C     .                   chi2 = chi2 + ((93.0d0 - mh)*120+ 200)*fac
C      if (mh.lt.92.0d0.and.mh.ge.91.0d0) 
C     .                   chi2 = chi2 + ((mh - 92.0d0)*30 + 320)*fac
C      if (mh.lt.91.0d0.and.mh.ge.89.5d0) 
C     .                   chi2 = chi2 + ((91.0d0 - mh)*54 + 290)*fac
C      if (mh.lt.89.5d0)  chi2 = chi2 + ((89.5d0 - mh)*58 + 371)*fac

C  summer 1998 Higgs exclusion (LEP):

C      fac = 2*dlog(10.d0)/100
C      if (mh.lt.97.0d0.and.mh.ge.92.0d0) 
C     .                   chi2 = chi2 + (97 - mh)*10*fac
C      if (mh.lt.92.0d0.and.mh.ge.89.6d0) 
C     .                   chi2 = chi2 + ((92     - mh)*45 +  50)*fac
C      if (mh.lt.89.6d0.and.mh.ge.88.4d0) 
C     .                   chi2 = chi2 + ((89.6d0 - mh)*90 + 158)*fac
C      if (mh.lt.88.4d0.and.mh.ge.86.d0) 
C     .                   chi2 = chi2 + ((88.4d0 - mh)*60 + 266)*fac
C      if (mh.lt.86.0d0.and.mh.ge.84.4d0) 
C     .                   chi2 = chi2 + ((86     - mh)*85 + 410)*fac
C      if (mh.lt.84.4d0)  chi2 = chi2 + ((84.4d0 - mh)*35 + 546)*fac

C  winter 1998 Higgs exclusion (ALEPH):

C       fac = 2*dlog(10.d0)/100
C       if (mh.lt.96.d0.and.mh.ge.92d0) chi2 = chi2 + (96 - mh)*6*fac
C       if (mh.lt.92.d0.and.mh.ge.88d0) 
C     .      chi2 = chi2 + ((92 - mh)*27 +  24)*fac
C       if (mh.lt.88.d0.and.mh.ge.84d0) 
C     .      chi2 = chi2 + ((88 - mh)*15 + 132)*fac
C       if (mh.lt.84.d0.and.mh.ge.83d0) chi2 = chi2 +        192 *fac
C       if (mh.lt.83.d0.and.mh.ge.78d0) 
C     .      chi2 = chi2 + ((83 - mh)*26 + 192)*fac
C       if (mh.lt.78.d0)       chi2 = chi2 + ((78 - mh)*12 + 322)*fac
       
C  summer 1997 Higgs exclusion (LEP):

C      fac = 2*dlog(10.d0)
C       if (mh.lt.83.d0.and.mh.ge.74d0) chi2 = chi2 + (83 - mh)/9*2*fac
C       if (mh.lt.74.d0.and.mh.ge.73d0) chi2 = chi2 +             2*fac
C       if (mh.lt.73.d0)         chi2 = chi2 +  ((73 - mh)/7*5 + 2)*fac
