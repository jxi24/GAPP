      subroutine values

      implicit none
      integer i,j,k,f,l,numrun
      double precision gammaz(0:11),sigmah,R(0:9),sin2te(0:9),alr(0:10)
      double precision afb(0:10),dkapse,QwE158,QwJLab,C1u,C1d,C2u,C2d
      double precision alphas,QwH,QwCs,QwTl,QwFr,QFr207,QFr211,Qfr228,QB10
      double precision mtrun,mbrun,mcrun,msrun,mtpole,mcpole,mbpole
      double precision CCFR,NuTeV,CDHS1,CDHS2,CDHS3,gvnue,ganue
      double precision epsu_L,epsd_L,epsu_R,epsd_R,sgcenu,gammaw(7)
      double precision z(0:13,80),rho(12,12),emat(12,12)
      double precision eplus(20),eminus(20),eparab(20),globcc(20)
      double precision alfame,alfmmu,alfaud,alfams,amcbar,alfamc,almtau
      double precision alfajpsi,alfamb,alfaups,alfamw,alfamt,alfax
      double precision dgamma,taulif,pterr2,amu,shat,dummy
      complex*16 kappa(0:9),rhobar(0:9)
      include 'common.f'

      call mnemat(emat,12)

      do 3 l = 0, 12
         if (l.gt.0) 
     .        call mnerrs(l,eplus(l),eminus(l),eparab(l),globcc(l))
         if (l.eq. 1) mz     = mz     + eparab(1)
         if (l.eq. 2) mt     = mt     + eparab(2)
         if (l.eq. 3) mb     = mb     + eparab(3)
         if (l.eq. 4) mc     = mc     + eparab(4)
         if (l.eq. 5) alfas0 = alfas0 + eparab(5)
         if (l.eq. 6) dahad3 = dahad3 + eparab(6)
         if (l.eq. 7) mh     =  mh*dexp(eparab(7))
         if (l.eq. 8) Tpar   = Tpar   + eparab(8)
         if (l.eq. 9) Spar   = Spar   + eparab(9)
         if (l.eq.10) Upar   = Upar   + eparab(10)
         if (l.eq.11) Brho   = Brho   + eparab(11)
         if (l.eq.12) then
            Bkappa = Bkappa + eparab(12)
            do 2 i = 1, 12
               do 1 j = 1, 12
                  rho(i,j) = 0.d0
                  if (eparab(i).ne.0.d0.and.eparab(j).ne.0.d0)
     .                 rho(i,j) = emat(i,j)/eparab(i)/eparab(j)
 1             continue
 2          continue
         endif

         mz2 = mz**2
         mt2 = mt**2
         mh2 = mh**2

      rathz2 = mh2/mz2
      rattz2 = mt2/mz2
      ratth2 = mt2/mh2      

         call sin2thetaw
         call z0pole(gammaz,sigmah,R,sin2te,alr,afb,0,9)
         call wwprod(gammaw)
         call rho0
         call sin2theta0
         call apv(QwH , 1,  1,C1u,C1d,C2u,C2d)
         call apv(QwCs,55,133,C1u,C1d,C2u,C2d)
         call apv(QwTl,81,205,C1u,C1d,C2u,C2d)
         call bsgamma(sgcenu)
c     call apv(QwFr,87,210,C1u,C1d,C2u,C2d)
c     call apv(QFr207,87,207,C1u,C1d,C2u,C2d)
c     call apv(QFr211,87,211,C1u,C1d,C2u,C2d)
c     call apv(QFr228,87,228,C1u,C1d,C2u,C2d)
c     call apv(QB10,5,10,C1u,C1d,C2u,C2d)
         call taulifetime(taulif,pterr2)
         call anomagmntmu(amu)
         call nuhnutev(NuTeV,epsu_L,epsd_L,epsu_R,epsd_R)
         call nuhccfr(CCFR)
         call nuhcdhs(CDHS1,CDHS2,CDHS3)
         call nue(0.d0,gvnue,ganue)
         call moller(QwE158,0.0260d0,0.599d0)
         call moller(QwJLab,0.0056d0,0.571d0)

         call rhof(rhobar,0,9)
         call kappaf(kappa,0,9)
         call polemasses(6,mtpole)
         call polemasses(5,mbpole)
         call polemasses(4,mcpole)

         call alfahat(me,dgamma,alfame)
         call alfahat(mmu,dgamma,alfmmu)
         call alfahat(0.176d0,dgamma,alfaud)
         call alfahat(0.305d0,dgamma,alfams)
         call alfahat(1.176d0,dgamma,amcbar)
         call alfahat(mc,dgamma,alfamc)
         call alfahat(3.096870d0,dgamma,alfajpsi)
         call alfahat(mtau,dgamma,almtau)
         call alfahat(mb,dgamma,alfamb)
         call alfahat(9.4603d0,dgamma,alfaups)
         call alfahat(mw,dgamma,alfamw)
         call alfahat(mt,dgamma,alfamt)
         call alfahat(0.194936d0,dgamma,alfax)

         if (l.eq.0) then
C     C     write(7,*) 'flagmr          = ', flagmr
C     C     write(7,*) 'flgblm          = ', flgblm
C     C     write(7,*) 'flgech          = ', flgech
C     C     write(7,*) 'f4lqcd          = ', f4lqcd
C     C     write(7,*) 'ffermi          = ', ffermi
C     C     write(7,*) 'fa2mt4          = ', fa2mt4
C     C     write(7,*) 'fa2mt2          = ', fa2mt2
C     C     write(7,*) 'fa2mt0          = ', fa2mt0
C     C     write(7,*) 'fla2im          = ', fla2im
C     C     write(7,*) 'fasmt2          = ', fasmt2
C     C     write(7,*) 'fasmt0          = ', fasmt0
C     C     write(7,*) 'fas2mt          = ', fas2mt
C     C     write(7,*) 'falas2          = ', falas2
C     C     write(7,*) 'flagmf          = ', flagmf
C     C     write(7,*) 'fbayes          = ', fbayes
C     C     write(7,*) 'fobliq          = ', fobliq
C     C     write(7,*) 'flagzp          = ', flagzp
C     C     write(7,*)
C     C     write(7,*) 'S               = ', Spar
C     C     write(7,*) 'T               = ', Tpar
C     C     write(7,*) 'U               = ', Upar
C     C     write(7,*)
C     C     write(7,*) 'M_H             = ', mh
C     C     write(7,*) 'M_Z             = ', mz
C     C     write(7,*) 'delta alp^(3)   = ', dahad3
C     C     write(7,*)
C     C     write(7,*) 'm_t(1 GeV)       = ', mtrun(1.d0)
C     C     write(7,*) 'm_t(M_Z)         = ', mtrun(mz)
C     C     write(7,*) 'm_t(m_t*M_Z)     = ', mtrun(dsqrt(mz*mt))
C     C     write(7,*) 'm_t(m_t)         = ', mt
C     C     write(7,*) 'm_t(m_t^pole)    = ', mtrun(mtpole)
C     C     write(7,*) 'm_t^pole         = ', mtpole
C     C     write(7,*)
C     C     write(7,*) 'm_b(1 GeV)       = ', mbrun(1.d0)
C     C     write(7,*) 'm_b(2 GeV)       = ', mbrun(2.d0)
C     C     write(7,*) 'm_b(m_b)         = ', mb
C     C     write(7,*) 'm_b(Upsilon1S/2) = ', mbrun(4.73015d0)
C     C     write(7,*) 'm_b(M_B^+-)      = ', mbrun(5.279d0)
C     C     write(7,*) 'm_b(M_B*)        = ', mbrun(5.325d0)
C     C     write(7,*) 'm_b(Upsilon_1S)  = ', mbrun(9.4603d0)
C     C     write(7,*) 'm_b(2 M_B^+-  )  = ', mbrun(10.558d0)
C     C     write(7,*) 'm_b(M_Z)         = ', mbrun(mz)
C     C     write(7,*) 'm_b(1 TeV)       = ', mbrun(1.d3)
C     C     write(7,*) 'm_b^pole         = ', mbpole
C     C     write(7,*)
C     C     write(7,*) 'm_c(1 GeV)       = ', mcrun(1.d0)
C     C     write(7,*) 'm_c(2 GeV)       = ', mcrun(2.d0)
C     C     write(7,*) 'm_c(m_c)         = ', mc
C     C     write(7,*) 'm_c(J/Psi/2)     = ', mcrun(1.548435d0)
C     C     write(7,*) 'm_c(M_D^0)       = ', mcrun(1.8645d0)
C     C     write(7,*) 'm_c(M_D^+-)      = ', mcrun(1.8693d0)
C     C     write(7,*) 'm_c(M_D_s)       = ', mcrun(1.9685d0)
C     C     write(7,*) 'm_c(M_D^+-*)     = ', mcrun(2.0100d0)
C     C     write(7,*) 'm_c(M_D_s*)      = ', mcrun(2.1124d0)
C     C     write(7,*) 'm_c(J/Psi)       = ', mcrun(3.096870d0)
C     C     write(7,*) 'm_c(2 M_D^0)     = ', mcrun(3.729d0)
C     C     write(7,*) 'm_c(m_b)         = ', mcrun(mb)
C     C     write(7,*) 'm_c(M_Z)         = ', mcrun(mz)
C     C     write(7,*) 'm_c(1 TeV)       = ', mcrun(1.d3)
C     C     write(7,*) 'm_c^pole         = ', mcpole
C     C     write(7,*)
C     C     write(7,*) 'm_s(1000 MeV)    = ', msrun(1.d0)
C     C     write(7,*) 'm_d(1000 MeV)    = ', md
C     C     write(7,*) 'm_u(1000 MeV)    = ', mu 
C     C     write(7,*)
C     C     write(7,*) 'm_s(mtau)        = ', msrun(mtau)
C     C     write(7,*) 'm_d(mtau)        = ', md/ms*msrun(mtau)
C     C     write(7,*) 'm_u(mtau)        = ', mu/ms*msrun(mtau)
C     C     write(7,*)
C     C     write(7,*) 'm_s(2000 MeV)    = ', msrun(2.d0)
C     C     write(7,*)
C     C     write(7,*) 'alfas0           = ', alfas0
C     C     write(7,*) 'alphas(m_t)      = ', alphas(mt)
C     C     write(7,*) 'alphas(M_Z)      = ', alphas(mz)
C     C     write(7,*) 'alphas(M_W)      = ', alphas(mw)
C     C     write(7,*) 'alphas(2 M_B^+-) = ', alphas(10.558d0)
C     C     write(7,*) 'alphas(Upsilon)  = ', alphas(9.46030d0)
C     C     write(7,*) 'alphas(m_b)      = ', alphas(mb)
C     C     write(7,*) 'alphas(2 M_D^0)  = ', alphas(3.729d0)
C     C     write(7,*) 'alphas(J/Psi)    = ', alphas(3.09687d0)
C     C     write(7,*) 'alphas(1800 MeV) = ', alphas(1.8d0)
C     C     write(7,*) 'alphas(m_tau)    = ', alphas(mtau)
C     C     write(7,*) 'alphas(m_c)      = ', alphas(mc)
C     C     write(7,*) 'alphas(1176 MeV) = ', alphas(1.176d0)
C     C     write(7,*) 'alphas(Phi)      = ', alphas(1.019456d0)
C     C     write(7,*) 'alphas(7500 MeV) = ', alphas(7.5d0)
C     C     write(7,*)
C alpha_s values relevant for BES data:
C
CC     C     write(7,*) '2 alphas(2.000) = ', 2*alphas(2.000d0)/pi1
CC     C     write(7,*) '2 alphas(2.200) = ', 2*alphas(2.200d0)/pi1
CC     C     write(7,*) '2 alphas(2.400) = ', 2*alphas(2.400d0)/pi1
CC     C     write(7,*) '2 alphas(2.500) = ', 2*alphas(2.500d0)/pi1
CC     C     write(7,*) '2 alphas(2.600) = ', 2*alphas(2.600d0)/pi1
CC     C     write(7,*) '2 alphas(2.700) = ', 2*alphas(2.700d0)/pi1
CC     C     write(7,*) '2 alphas(2.800) = ', 2*alphas(2.800d0)/pi1
CC     C     write(7,*) '2 alphas(2.900) = ', 2*alphas(2.900d0)/pi1
CC     C     write(7,*) '2 alphas(3.000) = ', 2*alphas(3.000d0)/pi1
CC     C     write(7,*) '2 alphas(3.700) = ', 2*alphas(3.700d0)/pi1
CC     C     write(7,*) '2 alphas(3.730) = ', 2*alphas(3.730d0)/pi1
CC     C     write(7,*) '2 alphas(3.750) = ', 2*alphas(3.750d0)/pi1
CC     C     write(7,*) '2 alphas(3.760) = ', 2*alphas(3.760d0)/pi1
CC     C     write(7,*) '2 alphas(3.764) = ', 2*alphas(3.764d0)/pi1
CC     C     write(7,*) '2 alphas(3.768) = ', 2*alphas(3.768d0)/pi1
CC     C     write(7,*) '2 alphas(3.770) = ', 2*alphas(3.770d0)/pi1
CC     C     write(7,*) '2 alphas(3.772) = ', 2*alphas(3.772d0)/pi1
CC     C     write(7,*) '2 alphas(3.776) = ', 2*alphas(3.776d0)/pi1
CC     C     write(7,*) '2 alphas(3.780) = ', 2*alphas(3.780d0)/pi1
CC     C     write(7,*) '2 alphas(3.790) = ', 2*alphas(3.790d0)/pi1
CC     C     write(7,*) '2 alphas(3.810) = ', 2*alphas(3.810d0)/pi1
CC     C     write(7,*) '2 alphas(3.850) = ', 2*alphas(3.850d0)/pi1
CC     C     write(7,*) '2 alphas(3.890) = ', 2*alphas(3.890d0)/pi1
CC     C     write(7,*) '2 alphas(3.930) = ', 2*alphas(3.930d0)/pi1
CC     C     write(7,*) '2 alphas(3.940) = ', 2*alphas(3.940d0)/pi1
CC     C     write(7,*) '2 alphas(3.950) = ', 2*alphas(3.950d0)/pi1
CC     C     write(7,*) '2 alphas(3.960) = ', 2*alphas(3.960d0)/pi1
CC     C     write(7,*) '2 alphas(3.970) = ', 2*alphas(3.970d0)/pi1
CC     C     write(7,*) '2 alphas(3.980) = ', 2*alphas(3.980d0)/pi1
CC     C     write(7,*) '2 alphas(3.990) = ', 2*alphas(3.990d0)/pi1
CC     C     write(7,*) '2 alphas(4.000) = ', 2*alphas(4.000d0)/pi1
CC     C     write(7,*) '2 alphas(4.010) = ', 2*alphas(4.010d0)/pi1
CC     C     write(7,*) '2 alphas(4.020) = ', 2*alphas(4.020d0)/pi1
CC     C     write(7,*) '2 alphas(4.027) = ', 2*alphas(4.027d0)/pi1
CC     C     write(7,*) '2 alphas(4.030) = ', 2*alphas(4.030d0)/pi1
CC     C     write(7,*) '2 alphas(4.033) = ', 2*alphas(4.033d0)/pi1
CC     C     write(7,*) '2 alphas(4.040) = ', 2*alphas(4.040d0)/pi1
CC     C     write(7,*) '2 alphas(4.050) = ', 2*alphas(4.050d0)/pi1
CC     C     write(7,*) '2 alphas(4.060) = ', 2*alphas(4.060d0)/pi1
CC     C     write(7,*) '2 alphas(4.070) = ', 2*alphas(4.070d0)/pi1
CC     C     write(7,*) '2 alphas(4.080) = ', 2*alphas(4.080d0)/pi1
CC     C     write(7,*) '2 alphas(4.090) = ', 2*alphas(4.090d0)/pi1
CC     C     write(7,*) '2 alphas(4.100) = ', 2*alphas(4.100d0)/pi1
CC     C     write(7,*) '2 alphas(4.110) = ', 2*alphas(4.110d0)/pi1
CC     C     write(7,*) '2 alphas(4.120) = ', 2*alphas(4.120d0)/pi1
CC     C     write(7,*) '2 alphas(4.130) = ', 2*alphas(4.130d0)/pi1
CC     C     write(7,*) '2 alphas(4.140) = ', 2*alphas(4.140d0)/pi1
CC     C     write(7,*) '2 alphas(4.150) = ', 2*alphas(4.150d0)/pi1
CC     C     write(7,*) '2 alphas(4.160) = ', 2*alphas(4.160d0)/pi1
CC     C     write(7,*) '2 alphas(4.170) = ', 2*alphas(4.170d0)/pi1
CC     C     write(7,*) '2 alphas(4.180) = ', 2*alphas(4.180d0)/pi1
CC     C     write(7,*) '2 alphas(4.190) = ', 2*alphas(4.190d0)/pi1
CC     C     write(7,*) '2 alphas(4.200) = ', 2*alphas(4.200d0)/pi1
CC     C     write(7,*) '2 alphas(4.210) = ', 2*alphas(4.210d0)/pi1
CC     C     write(7,*) '2 alphas(4.220) = ', 2*alphas(4.220d0)/pi1
CC     C     write(7,*) '2 alphas(4.230) = ', 2*alphas(4.230d0)/pi1
CC     C     write(7,*) '2 alphas(4.240) = ', 2*alphas(4.240d0)/pi1
CC     C     write(7,*) '2 alphas(4.245) = ', 2*alphas(4.245d0)/pi1
CC     C     write(7,*) '2 alphas(4.250) = ', 2*alphas(4.250d0)/pi1
CC     C     write(7,*) '2 alphas(4.255) = ', 2*alphas(4.255d0)/pi1
CC     C     write(7,*) '2 alphas(4.260) = ', 2*alphas(4.260d0)/pi1
CC     C     write(7,*) '2 alphas(4.265) = ', 2*alphas(4.265d0)/pi1
CC     C     write(7,*) '2 alphas(4.270) = ', 2*alphas(4.270d0)/pi1
CC     C     write(7,*) '2 alphas(4.280) = ', 2*alphas(4.280d0)/pi1
CC     C     write(7,*) '2 alphas(4.300) = ', 2*alphas(4.300d0)/pi1
CC     C     write(7,*) '2 alphas(4.320) = ', 2*alphas(4.320d0)/pi1
CC     C     write(7,*) '2 alphas(4.340) = ', 2*alphas(4.340d0)/pi1
CC     C     write(7,*) '2 alphas(4.350) = ', 2*alphas(4.350d0)/pi1
CC     C     write(7,*) '2 alphas(4.360) = ', 2*alphas(4.360d0)/pi1
CC     C     write(7,*) '2 alphas(4.380) = ', 2*alphas(4.380d0)/pi1
CC     C     write(7,*) '2 alphas(4.390) = ', 2*alphas(4.390d0)/pi1
CC     C     write(7,*) '2 alphas(4.400) = ', 2*alphas(4.400d0)/pi1
CC     C     write(7,*) '2 alphas(4.410) = ', 2*alphas(4.410d0)/pi1
CC     C     write(7,*) '2 alphas(4.420) = ', 2*alphas(4.420d0)/pi1
CC     C     write(7,*) '2 alphas(4.430) = ', 2*alphas(4.430d0)/pi1
CC     C     write(7,*) '2 alphas(4.440) = ', 2*alphas(4.440d0)/pi1
CC     C     write(7,*) '2 alphas(4.450) = ', 2*alphas(4.450d0)/pi1
CC     C     write(7,*) '2 alphas(4.460) = ', 2*alphas(4.460d0)/pi1
CC     C     write(7,*) '2 alphas(4.480) = ', 2*alphas(4.480d0)/pi1
CC     C     write(7,*) '2 alphas(4.500) = ', 2*alphas(4.500d0)/pi1
CC     C     write(7,*) '2 alphas(4.520) = ', 2*alphas(4.520d0)/pi1
CC     C     write(7,*) '2 alphas(4.540) = ', 2*alphas(4.540d0)/pi1
CC     C     write(7,*) '2 alphas(4.560) = ', 2*alphas(4.560d0)/pi1
CC     C     write(7,*) '2 alphas(4.600) = ', 2*alphas(4.600d0)/pi1
CC     C     write(7,*) '2 alphas(4.800) = ', 2*alphas(4.800d0)/pi1

cC     C     write(7,*) 'Gamma(Upsilon --> gamma g g)/Gamma(g g g):'

cC     C     write(7,*) 0.0275d0,' +- ',0.0016d0,'     ',
c     .           4*alfaups/alphas(9.4603d0)/5*
c     .           (1 - 4.86d0*alphas(9.4603d0)/pi1),' +- ',4*alfaups/
c     .           alphas(9.4603d0)/5*(4.86d0*alphas(9.4603d0)/pi1)**2

cC     C     write(7,*) 'Gamma(Upsilon --> gamma g g)/Gamma(l^+ l^-):'

cC     C     write(7,*) 0.894d0,' +- ',0.062d0,'     ',
c     .           8*(pi2/9 - 1)*alphas(9.4603d0)**2/pi1/alfaups*
c     .           (1 + 2.43d0*alphas(9.4603d0)/pi1)*0.69d0,' +- ',
c     .        8*(pi2/9 - 1)*alphas(9.4603d0)**4/pi1/alfaups*9/pi2*0.69d0

cC     C     write(7,*) 'Gamma(J/Psi --> gamma g g)/Gamma(g g g):'

cC     C     write(7,*) 0.10d0,' +- ',0.04d0,'     ',
c     .           16*alfajpsi/alphas(3.09687d0)/5*
c     .           (1 - 5.48d0*alphas(3.09687d0)/pi1),' +- ',16*alfajpsi/
c     .           alphas(3.09687d0)/5*(5.48d0*alphas(3.09687d0)/pi1)**2

cC     C     write(7,*) 'Gamma(J/Psi --> gamma g g)/Gamma(l^+ l^-):'

cC     C     write(7,*) 0.90d0,' +- ',0.38d0,'     ',
c     .           8*(pi2/9 - 1)*alphas(3.09687d0)**2/pi1/alfajpsi*
c     .           (1 + 3.66d0*alphas(3.09687d0)/pi1)*0.31d0,' +- ',
c     .           8*(pi2/9 -1)*alphas(3.09687d0)**4/pi1/alfajpsi
c     .           *3.66d0**2/pi2*0.31d0

cC     C     write(7,*)

cC     C     write(7,*) 10*(pi2 - 9)*alphas(9.4603d0)**3/9/pi1/alfaups**2
c     .          *(1 + 7.28d0*alphas(9.4603d0)/pi1)*0.69d0,
c     .           10*(pi2 - 9)*alphas(9.4603d0)**3/9/pi1/alfaups**2
c     .          *(7.28d0*alphas(9.4603d0)/pi1)**2*0.69d0

cC     C     write(7,*) 5*(pi2 - 9)*alphas(3.09687d0)**3/18/pi1/
c     .           alfajpsi**2*(1 + 9.13d0*alphas(3.09687d0)/pi1)*0.31d0,
c     .           5*(pi2 - 9)*alphas(3.09687d0)**3/18/pi1/
c     .           alfajpsi**2*(9.13d0*alphas(3.09687d0)/pi1)**2*0.31d0
 
cC     C     write(7,*) 4*(alphas(9.4603d0)/alphas(3.09687d0))**3*
c     .           (alfamc/alfaups)**2*(1 + 7.28d0*alphas(9.4603d0)/pi1 
c     .           - 9.13d0*alphas(3.09687d0)/pi1)/0.45d0

cC     C     write(7,*) 4*(alphas(9.4603d0)/alphas(3.09687d0))**3*
c     .           (alfamc/alfaups)**2*(1 + 7.28d0*alphas(9.4603d0)/pi1)/ 
c     .           (1 + 9.13d0*alphas(3.09687d0)/pi1)/0.45d0

cC     C     write(7,*)

cC     C     write(7,*) 16*alfajpsi/alphas(3.09687d0)/5*
c     .           (1 - 5.48d0*alphas(3.09687d0)/pi1),16*alfajpsi/
c     .           alphas(3.09687d0)/5*(5.48d0*alphas(3.09687d0)/pi1)**2

cC     C     write(7,*) alfaups/alfajpsi/alphas(9.4603d0)/4*
c     .           alphas(3.09687d0)*(1 - 4.86d0*alphas(9.4603d0)/pi1)/
c     .           (1 - 5.48d0*alphas(3.09687d0)/pi1)

cC     C     write(7,*)

cC     C     write(7,*) 10*(pi2 - 9)*alphas(mb)**3/9/pi1/alfamb**2*
c     .           (1 + 0.47d0*alphas(mb)/pi1)*0.69d0

cC     C     write(7,*) 5*(pi2 - 9)*alphas(mc)**3/18/pi1/alfamc**2*
c     .           (1 + 1.63d0*alphas(mc)/pi1)*0.31d0

cC     C     write(7,*) 4*(alphas(mb)/alphas(mc))**3*
c     .           (alfamc/alfaups)**2*(1 - 0.73d0*alphas(mb)/pi1 
c     .           - 0.43d0*alphas(mc)/pi1)/0.45d0

cC     C     write(7,*) 4*(alphas(mb)/alphas(mc))**3*
c     .           (alfamc/alfaups)**2*(1 + 0.47d0*alphas(mb)/pi1)/ 
c     .           (1 + 1.63d0*alphas(mc)/pi1)/0.45d0

cC     C     write(7,*)

cC     C     write(7,*) 4*alfaups/alphas(mb)/5*
c     .           (1 - 2.16d0*alphas(mb)/pi1)

cC     C     write(7,*) 16*alfajpsi/alphas(mc)/5*
c     .           (1 - 2.55d0*alphas(mc)/pi1)

cC     C     write(7,*) alfaups/alfajpsi/alphas(mb)*alphas(mc)/4*
c     .           (1 - 2.16d0*alphas(mb)/pi1)/(1 - 2.55d0*alphas(mc)/pi1)

C  For the extraction of Lambda_QCD (3-flavor; 4-loop):
C
C      do 99 i = 1,900
C   C     write(7,*) alphas(0.300d0 + i*0.001d0),i
C 99   continue

C     C     write(7,*) '1/alpha           = ', 1/alpha
C     C     write(7,*) '1/alpha(m_e)      = ', 1/alfame
C     C     write(7,*) '1/alpha(m_mu)     = ', 1/alfmmu
C     C     write(7,*) '1/alpha( 176 MeV) = ', 1/alfaud
C     C     write(7,*) '1/alpha( 305 MeV) = ', 1/alfams
C     C     write(7,*) '1/alpha(1176 MeV) = ', 1/amcbar
C     C     write(7,*) '1/alpha(m_c)      = ', 1/alfamc
C     C     write(7,*) '1/alpha(m_tau)    = ', 1/almtau
C     C     write(7,*) '1/alpha(J/Psi)    = ', 1/alfajpsi
C     C     write(7,*) '1/alpha(m_b)      = ', 1/alfamb
C     C     write(7,*) '1/alpha(Upsilon)  = ', 1/alfaups
C     C     write(7,*) '1/alpha(M_W)      = ', 1/alfamw
C     C     write(7,*) '1/alpha(M_Z)      = ', 1/alphat
C     C     write(7,*) '1/alpha(x)        = ', 1/alfax
C     C     write(7,*) '1/alpha(m_t)      = ', 1/alfamt
C     C     write(7,*) 'Delta alphat      = ', 1 - alpha/alphat
C     C     write(7,*)
C     C     write(7,*) 'shat(M_GUT)        = ', shat(1.d13)
C     C     write(7,*) 'shat(1960 GeV)     = ', shat(1.96d3)
C     C     write(7,*) 'shat(mt)           = ', shat(mt)
C     C     write(7,*) 'shat(x)            = ', shat(mt - 1.d-12)
C     C     write(7,*) 'shat(mz)           = ', shat(mz)
C     C     write(7,*) 'shat(mw)           = ', shat(mw)
C     C     write(7,*) 'shat(mb)           = ', shat(mb)
C     C     write(7,*) 'shat(20 GeV^2)     = ', shat(4.472d0)
C     C     write(7,*) 'shat(6.25 GeV^2)   = ', shat(2.5d0)
C     C     write(7,*) 'shat(mtau)         = ', shat(mtau)
C     C     write(7,*) 'shat(mc)           = ', shat(mc)
C     C     write(7,*) 'shat(1176 MeV)     = ', shat(1.176d0)
C     C     write(7,*) 'shat( 305 MeV)     = ', shat(0.305d0)
C     C     write(7,*) 'shat( 176 MeV)     = ', shat(0.176d0)
C     C     write(7,*) 'shat(0.026 GeV^2)  = ', shat(0.161d0)
C     C     write(7,*) 'shat(mmu)          = ', shat(mmu)
C     C     write(7,*) 'shat(0.0056 GeV^2) = ', shat(0.075d0)
C     C     write(7,*) 'shat(me)           = ', shat(me)
C     C     write(7,*) 'shat(   0 MeV)     = ', shat(0.d0)
C     C     write(7,*) 'sin^2th(M_Z)       = ', sinhat
C     C     write(7,*) 'sin^2th_W          = ', 1 - mw2/mz2
C     C     write(7,*) 'M_W                = ', mw
C     C     write(7,*)
C     C     write(7,*) 'GammaW_enu      = ', gammaw(1)
C     C     write(7,*) 'GammaW_ud/V_ud  = ', gammaw(4)
C     C     write(7,*) 'GammaW_had      = ', gammaw(6)
C     C     write(7,*) 'Gamma_W         = ', gammaw(7)
C     C     write(7,*)
C     C     write(7,*) 'GammaZ_inv      = ', gammaz(0)
C     C     write(7,*) 'GammaZ_e        = ', gammaz(1)
C     C     write(7,*) 'GammaZ_mu       = ', gammaz(2)
C     C     write(7,*) 'GammaZ_tau      = ', gammaz(3)
C     C     write(7,*) 'GammaZ_u        = ', gammaz(4)
C     C     write(7,*) 'GammaZ_c        = ', gammaz(5)
C     C     write(7,*) 'GammaZ_d        = ', gammaz(7)
C     C     write(7,*) 'GammaZ_s        = ', gammaz(8)
C     C     write(7,*) 'GammaZ_b        = ', gammaz(9)
C     C     write(7,*) 'GammaZ_had      = ', gammaz(10)
C     C     write(7,*) 'Gamma_Z         = ', gammaz(11)
C     C     write(7,*)
C     C     write(7,*) 'sin^eff_e       = ', sin2te(1)
C     C     write(7,*) 'sin^eff_c       = ', sin2te(5)
C     C     write(7,*) 'sin^eff_b       = ', sin2te(9)
C     C     write(7,*) 'R_e             = ', R(1)
C     C     write(7,*) 'R_mu            = ', R(2)
C     C     write(7,*) 'R_tau           = ', R(3)
C     C     write(7,*) 'R_u             = ', R(4)
C     C     write(7,*) 'R_d             = ', R(7)
C     C     write(7,*) 'R_s             = ', R(8)
C     C     write(7,*) 'R_c             = ', R(5)
C     C     write(7,*) 'R_b             = ', R(9)
C     C     write(7,*) 'sigma^0_h       = ', sigmah
C     C     write(7,*) 'A_FB^e          = ', afb(1)
C     C     write(7,*) 'A_FB^c          = ', afb(5)
C     C     write(7,*) 'A_FB^s          = ', afb(8)
C     C     write(7,*) 'A_FB^b          = ', afb(9)
C     C     write(7,*) 'A_FB^h          = ', afb(10)
C     C     write(7,*) 'A_FB^(u,d)      = ', afb(4)
C     C     write(7,*) 'A_LR            = ', alr(1)
C     C     write(7,*) 'A_c             = ', alr(5)
C     C     write(7,*) 'A_b             = ', alr(9)
C     C     write(7,*) 'A_s             = ', alr(8)
C     C     write(7,*)
C     C     write(7,*) 'rho_NC          = ', rhonc
C     C     write(7,*) 'sin^2th(0)      = ', shat(0.d0)
C     C     write(7,*) 'C_1u            = ', C1u
C     C     write(7,*) 'C_1d            = ', C1d
C     C     write(7,*) 'C_2u            = ', C2u
C     C     write(7,*) 'C_2d            = ', C2d
C     C     write(7,*) 'C_2u - 1/2 C_2d = ', C2u - 1/2.d0*C2d
C     C     write(7,*) '(2*C1u - C1d) + 0.84d0*(2*C2u - C2d) = ', 
C    .           (2*C1u - C1d) + 0.84d0*(2*C2u - C2d)
C     C     write(7,*) 'C_1u + C_1d     = ', C1u + C1d
C     C     write(7,*) 'C_1u - C_1d     = ', C1u - C1d
C     C     write(7,*) 'C_2u + C_2d     = ', C2u + C2d
C     C     write(7,*) 'C_2u - C_2d     = ', C2u - C2d
C     C     write(7,*) 'Q_W (e) [E158]  = ', QwE158
C     C     write(7,*) 'Q_W (e) [JLab]  = ', QwJLab
C     C     write(7,*) 'Q_W (H)         = ', QwH
C     C     write(7,*) 'Q_W (Cs)        = ', QwCs
C     C     write(7,*) 'Q_W (Tl)        = ', QwTl
c     write(7,*) 'Q_W (Fr)      = ', QwFr
c     write(7,*) 'Q_W (Fr207)   = ', QFr207
c     write(7,*) 'Q_W (Fr211)   = ', QFr211
c     write(7,*) 'Q_W (Fr228)   = ', QFr228
c     write(7,*) 'Q_W rat1      = ', QFr207/QFr211
c     write(7,*) 'Q_W rat2      = ', QFr207/QFr228
c     write(7,*) 'Q_W (B10)     = ', QB10
c     write(7,*) 'A_Fr(211-207) = ', (QFr211 - QFr207)/(QFr211 + QFr207)
c     write(7,*) 'A_Fr(228-207) = ', (QFr228 - QFr207)/(QFr228 + QFr207)
C     C     write(7,*) 'R_bsgamma       = ', sgcenu
C     C     write(7,*) 'taulifetime     = ', taulif
C     C     write(7,*)
C     C     write(7,*) 'kappa_e         = ', kappa(1)
C     C     write(7,*) 'kappa_u         = ', kappa(4)
C     C     write(7,*) 'kappa_c         = ', kappa(5)
C     C     write(7,*) 'kappa_d         = ', kappa(7)
C     C     write(7,*) 'kappa_s         = ', kappa(8)
C     C     write(7,*) 'kappa_b         = ', kappa(9)
C     C     write(7,*)
C     C     write(7,*) 'rho_e           = ', rhobar(1)
C     C     write(7,*) 'rho_u           = ', rhobar(4)
C     C     write(7,*) 'rho_c           = ', rhobar(5)
C     C     write(7,*) 'rho_d           = ', rhobar(7)
C     C     write(7,*) 'rho_s           = ', rhobar(8)
C     C     write(7,*) 'rho_b           = ', rhobar(9)
C     C     write(7,*)
C     C     write(7,*) 'kappa_CCFR      = ', CCFR
C     C     write(7,*) 'R^(-)_NuTeV     = ', NuTeV
C     C     write(7,*) 'R_nu_CDHS       = ', CDHS1
C     C     write(7,*) 'R_nubar_CDHS1   = ', CDHS2
C     C     write(7,*) 'R_nubar_CDHS2   = ', CDHS3
C     C     write(7,*)
C     C     write(7,*) 'g_V^nue         = ', gvnue
C     C     write(7,*) 'g_A^nue         = ', ganue
C     C     write(7,*)
C     C     write(7,*) 'epsilon_L(u)    = ', epsu_L
C     C     write(7,*) 'epsilon_L(d)    = ', epsd_L
C     C     write(7,*) 'epsilon_R(u)    = ', epsu_R
C     C     write(7,*) 'epsilon_R(d)    = ', epsd_R
C     C     write(7,*)
C     C     write(7,*) 'g_L^2           = ', epsu_L**2 + epsd_L**2
C     C     write(7,*) 'g_R^2           = ', epsu_R**2 + epsd_R**2
C     C     write(7,*) 'theta_L         = ', datan(epsu_L/epsd_L) +  pi1
C     C     write(7,*) 'theta_R         = ', datan(epsu_R/epsd_R) +2*pi1
C     C     write(7,*)
C     C     write(7,*) 'del_r           = ', delr
C     C     write(7,*) 'del_r_W         = ', delrwh
C     C     write(7,*) 'del_r_Z         = ', delrzh
C     C     write(7,*) 'rho_(MS-bar)    = ', rhohat
C     C     write(7,*)
         endif

         z(l, 1) = mtpole
         z(l, 2) = mh
         z(l, 3) = mw
         z(l, 4) = gammaw(7)
         z(l, 5) = gammaw(1)
         z(l, 6) = gammaw(4)
         z(l, 7) = mz
         z(l, 8) = gammaz(11)
         z(l, 9) = gammaz(10)
         z(l,10) = gammaz(0)
         z(l,11) = gammaz(1)
         z(l,12) = gammaz(4)
         z(l,13) = gammaz(7)
         z(l,14) = gammaz(9)
         z(l,15) = sigmah
         z(l,16) = R(1)
         z(l,17) = R(2)
         z(l,18) = R(3)
         z(l,19) = R(9)
         z(l,20) = R(5)
         z(l,21) = R(8)/(R(4)+R(7)+R(8))
         z(l,22) = afb(10)
         z(l,23) = afb(1)
         z(l,24) = afb(9)
         z(l,25) = afb(5)
         z(l,26) = afb(8)
         z(l,27) = alr(1)
         z(l,28) = alr(9)
         z(l,29) = alr(5)
         z(l,30) = alr(8)
         z(l,31) = 1 - mw2/mz2
         z(l,32) = NuTeV
         z(l,33) = CCFR
         z(l,34) = CDHS1
         z(l,35) = CDHS2
         z(l,36) = CDHS3
         z(l,37) = gvnue
         z(l,38) = ganue
         z(l,39) = QwE158
         z(l,40) = QwH
         z(l,41) = QwCs
         z(l,42) = QwTl
         z(l,43) = dlog(sgcenu)
         z(l,44) = delr
         z(l,45) = delrwh
         z(l,46) = epsu_L
         z(l,47) = epsd_L
         z(l,48) = epsu_R
         z(l,49) = epsd_R
         z(l,50) = epsu_L**2 + epsd_L**2
         z(l,51) = epsu_R**2 + epsd_R**2
         z(l,52) = datan(epsu_L/epsd_L) +   pi1
         z(l,53) = datan(epsu_R/epsd_R) + 2*pi1
         z(l,54) = C1u
         z(l,55) = C1d
         z(l,56) = C2u - 1/2.d0*C2d
         z(l,57) = C1u + C1d
         z(l,58) = C1u - C1d
         z(l,59) = C2u + C2d
         z(l,60) = C2u - C2d
         z(l,61) = 1/alphat
         z(l,62) = 1/almtau
         z(l,63) = sinhat
         z(l,64) = sin2te(1)
         z(l,65) = shat(0.d0)
         z(l,66) = taulif
         z(l,67) = alphas(mtau)
         z(l,68) = amu
         z(l,69) = shat(1.96d3)
         z(l,70) = 1/(1 - alphat*Tpar)
         z(l,71) = alphat*Spar/4/sinhat
         z(l,72) = alphat*Tpar
         z(l,73) =-alphat*Upar/4/sinhat
         z(l,74) = dsqrt(mz02)
         z(l,75) = sinth
         z(l,76) = ratg21
         z(l,77) = mzp
         z(l,78) = shat(mw)
         z(l,79) = mbrun(mb)
         z(l,80) = mt

         if (l.eq. 1) mz     = mz     - eparab(1)
         if (l.eq. 2) mt     = mt     - eparab(2)
         if (l.eq. 3) mb     = mb     - eparab(3)
         if (l.eq. 4) mc     = mc     - eparab(4)
         if (l.eq. 5) alfas0 = alfas0 - eparab(5)
         if (l.eq. 6) dahad3 = dahad3 - eparab(6)
         if (l.eq. 7) mh     =  mh/dexp(eparab(7))
         if (l.eq. 8) Tpar   = Tpar   - eparab(8)
         if (l.eq. 9) Spar   = Spar   - eparab(9)
         if (l.eq.10) Upar   = Upar   - eparab(10)
         if (l.eq.11) Brho   = Brho   - eparab(11)
         if (l.eq.12) Bkappa = Bkappa - eparab(12)
c         if (l.eq.13) mzp    = mzp/dexp(eparab(13))

 3    continue

      do 30 i = 1,80
         z(13,i) = 0.d0
         do 20 j = 1,12
            do 10 k = 1,12
               z(13,i)=z(13,i)+(z(j,i)-z(0,i))*rho(j,k)*(z(k,i)-z(0,i))
 10         continue

C   Error calculation with M_Z, M_H, m_t, and/or alpha(M_Z) fixed:
C   an inconvenience in MINUIT requires to neglect the error matrix;
C
C             z(13,i) = z(13,i)+(z(j,i)-z(0,i))*(z(j,i)-z(0,i))
C
C   alternatively, relabel fit parameters!

 20      continue
         z(13,i) = dsqrt(z(13,i))
 30   continue
      
C     write(7,*) '1 sigma errors:'
C     write(7,*) 
C     write(7,100) 'm_t (pole):      ', z(0, 1), z(13, 1) 
C     write(7,100) 'm_t (MS-bar):    ',      mt, eplus(2), eminus(2)
C     write(7,100) 'M_H (parabolic): ', z(0, 2), z(13, 2)
C     write(7,100) 'M_H (exact):     ', z(0, 2), z( 0, 2)*
C   .            (dexp(eplus(7)) - 1), z(0, 2)*(dexp(eminus(7)) - 1)
C     write(7,100) 'M_W:             ', z(0, 3), z(13, 3)
C     write(7,100) 'Gamma_W:         ', z(0, 4), z(13, 4)
C     write(7,100) 'GammaW_enu:      ', z(0, 5), z(13, 5)
C     write(7,100) 'GammaW_ud/V_ud:  ', z(0, 6), z(13, 6)
C     write(7,100) 'M_Z:             ', z(0, 7), z(13, 7)
C     write(7,100) 'Gamma_Z:         ', z(0, 8), z(13, 8)
C     write(7,100) 'GammaZ_had:      ', z(0, 9), z(13, 9)
C     write(7,100) 'GammaZ_inv:      ', z(0,10), z(13,10)
C     write(7,100) 'GammaZ_e:        ', z(0,11), z(13,11)
C     write(7,100) 'GammaZ_u:        ', z(0,12), z(13,12)
C     write(7,100) 'GammaZ_d:        ', z(0,13), z(13,13)
C     write(7,100) 'GammaZ_b:        ', z(0,14), z(13,14)
C     write(7,100) 'sigma_had:       ', z(0,15), z(13,15)
C     write(7,100) 'R_e:             ', z(0,16), z(13,16)
C     write(7,100) 'R_mu:            ', z(0,17), z(13,17)
C     write(7,100) 'R_tau:           ', z(0,18), z(13,18)
C     write(7,100) 'R_b:             ', z(0,19), z(13,19)
C     write(7,100) 'R_c:             ', z(0,20), z(13,20)
C     write(7,100) 'R_s/R_(u+d+s):   ', z(0,21), z(13,21)
C     write(7,100) 'A_FB(hadrons):   ', z(0,22), z(13,22)
C     write(7,100) 'A_FB(e):         ', z(0,23), z(13,23)
C     write(7,100) 'A_FB(b):         ', z(0,24), z(13,24)
C     write(7,100) 'A_FB(c):         ', z(0,25), z(13,25)
C     write(7,100) 'A_FB(s):         ', z(0,26), z(13,26)
C     write(7,100) 'A_LR:            ', z(0,27), z(13,27)
C     write(7,100) 'A_b:             ', z(0,28), z(13,28)
C     write(7,100) 'A_c:             ', z(0,29), z(13,29)
C     write(7,100) 'A_s:             ', z(0,30), z(13,30)
C     write(7,100) 'sin^2th_W:       ', z(0,31), z(13,31)
C     write(7,100) 'R^-:             ', z(0,32), z(13,32)
C     write(7,100) 'R^nu (CCFR):     ', z(0,33), z(13,33)
C     write(7,100) 'R^nu (CDHS):     ', z(0,34), z(13,34)
C     write(7,100) 'R^nubar (CDHS84):', z(0,35), z(13,35)
C     write(7,100) 'R^nubar (CDHS79):', z(0,36), z(13,36)
C     write(7,100) 'g_V^nue:         ', z(0,37), z(13,37)
C     write(7,100) 'g_A^nue:         ', z(0,38), z(13,38)
C     write(7,100) 'Q_W (e):         ', z(0,39), z(13,39)
C     write(7,100) 'Q_W (p):         ', z(0,40), z(13,40)
C     write(7,100) 'Q_W (Cs):        ', z(0,41), z(13,41)
C     write(7,100) 'Q_W (Tl):        ', z(0,42), z(13,42)
C     write(7,100) 'R_bsgamma:       ', z(0,43), z(13,43)
C     write(7,100) 'del_r:           ', z(0,44), z(13,44)
C     write(7,100) 'del_r_W:         ', z(0,45), z(13,45)
C     write(7,100) 'epsilon_L(u):    ', z(0,46), z(13,46)
C     write(7,100) 'epsilon_L(d):    ', z(0,47), z(13,47)
C     write(7,100) 'epsilon_R(u):    ', z(0,48), z(13,48)
C     write(7,100) 'epsilon_R(d):    ', z(0,49), z(13,49)
C     write(7,100) 'g_L^2:           ', z(0,50), z(13,50)
C     write(7,100) 'g_R^2:           ', z(0,51), z(13,51)
C     write(7,100) 'theta_L:         ', z(0,52), z(13,52)
C     write(7,100) 'theta_R:         ', z(0,53), z(13,53)
C     write(7,100) 'C_1u:            ', z(0,54), z(13,54)
C     write(7,100) 'C_1d:            ', z(0,55), z(13,55)
C     write(7,100) 'C_2u - C_2d/2:   ', z(0,56), z(13,56)
C     write(7,100) 'C_1u + C_1d:     ', z(0,57), z(13,57)
C     write(7,100) 'C_1u - C_1d:     ', z(0,58), z(13,58)
C     write(7,100) 'C_2u + C_2d:     ', z(0,59), z(13,59)
C     write(7,100) 'C_2u - C_2d:     ', z(0,60), z(13,60)
C     write(7,100) '1/alpha MS-bar:  ', z(0,61), z(13,61)
C     write(7,100) '1/alpha (m_tau): ', z(0,62), z(13,62)
C     write(7,100) 'sin^2 MS-bar:    ', z(0,63), z(13,63)
C     write(7,100) 'sin^2 effective: ', z(0,64), z(13,64)
C     write(7,100) 'sin^2 (0):       ', z(0,65), z(13,65)
C     write(7,100) 'taulifetime:     ', z(0,66), z(13,66)
C     write(7,100) 'alpha_s(m_tau):  ', z(0,67), z(13,67)
C     write(7,100) '(g-2-alfa/pi)/2: ', z(0,68), z(13,68)
C     write(7,100) 'sin^2 (1.96 TeV):', z(0,69), z(13,69)
c     write(7,100) 'rho_0 parameter: ', z(0,70), z(13,70)
C     write(7,100) 'rho_0 parameter: ', z(0,70), eparab(8)/z(0,61)
C     write(7,100) 'epsilon_3:       ', z(0,71), z(13,71)
C     write(7,100) 'epsilon_1:       ', z(0,72), z(13,72)
C     write(7,100) 'epsilon_2:       ', z(0,73), z(13,73)
C     write(7,100) 'M_Z_0:           ', z(0,74), z(13,74)
      call mnerrs(17,eplus(17),eminus(17),eparab(17),globcc(17))
C     write(7,100) 'sin(theta):      ', z(0,75), 
C    .     z(0,75) + eminus(17), z(0,75) + eplus(17)
      call mnerrs(16,eplus(16),eminus(16),eparab(16),globcc(16))
C     write(7,100) 'g_2/g_1:         ', z(0,76), z(13,76)
C     write(7,100) 'M_Z` (exact):    ', z(0,77),
C    .     z(0,77)*dexp(eminus(16)), z(0,77)*dexp(eplus(16))
C     write(7,100) 'M_Z` (exact):    ', mzp,
C    .     1/(1/mzp + 1.d-3*eminus(16)), 1/(1/mzp + 1.d-3*eplus(16))
C     write(7,100) 'placeholder:     ', z(0,78), z(13,78)
C     write(7,100) 'placeholder:     ', z(0,79), z(13,79)
C     write(7,100) 'placeholder:     ', z(0,80), z(13,80)
      
C  for rho_0 parameter:

      call mnerrs(8,eplus(8),eminus(8),eparab(8),globcc(8))
C     write(7,*) z(0,70),   eplus(8)/z(0,61), eminus(8)/z(0,61)
C     write(7,*) z(0,70) + eminus(8)/z(0,61), '< rho_0 <',
C    .           z(0,70) +  eplus(8)/z(0,61)

 100  format(a20,f10.5,f14.7,f14.7,f10.3)
      
      stop
      end
