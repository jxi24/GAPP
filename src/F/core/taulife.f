      subroutine taulifetime(taulif,therr2)

C  alpha^2 corrections to muon decays have been calculated by R. Stuart and 
C  T. van Ritbergen in hep-ph/9808283, 9802341, and 9904240.  The analogous 
C  effects for tau decays are about 3.6*10^-5 and can be obtained by means 
C  of a dispersion relation. We neglect effects of this order.
C
C  Phase space correction for muon final states including O(alpha) 
C  contribution extracted from Y. Nir, PLB 221 (1989) 184.
C
C  One-loop electroweak corrections from
C  W.J. Marciano and A. Sirlin, PRL 61, 1815 (1988), and
C  E. Braaten and C.S. Li, PRD 42, 3888 (1990).
C
C  Massless QCD corrections with resummation according to
C  F. de LeDiberder and A. Pich, PLB 286, 147 (1992).
C
C  Charm quark mass corrections in an expansion in (mtau/mc)**2 using
C  S.A. Larin, T. van Ritbergen, and J.A.M. Vermaseren, NPB 438, 278 (1995). 

      implicit none
      double precision alphas,as,astau,asb,asw,etat,etab,etaw,cmatch,c3m
      double precision arun,az,aw,ab,atau,dgamma
      double precision g0,ge,gmu,gud,gtot,vud,hbar,rmt2,taulif,therr2
      double precision sew,delqcd,pterr,remerr,rustot
      double precision cdim4,cdim6,cdim8,mpipm,mkpm,fpipm,fkpm
      include 'common.f'

      cmatch = 11/72.d0
        hbar = 6.58211889d-25
c        vud = 0.97451d0   ! +- 0.00037 using s_ij and delta_13 from PDG 2004
c        vud = 0.97418d0   ! +- 0.00027 using 0801.1821
         vud = 0.97408d0   ! +- 0.00026 using 0807.0650
       mpipm = 0.13957018d0
        mkpm = 0.493677d0
       fpipm = 0.0924d0
        fkpm = dsqrt(1.1d0)*fpipm ! +- sqrt(0.6)*fpipm
      rustot = 0.0285d0
C     remerr = 0.9d0 
      remerr = 1.1d0   ! needs update, also alpha_s from Z decays

C Non-perturbative contribution according to 0807.0650:

      cdim4  =   0.009d0      ! +- 0.007
      cdim6  = - 0.000254d0/2 ! +- 0.00030/2 (exp.) +- 0.000114/2 (th.) = 0.00059
      cdim8  =   0.000236d0/2 ! +- 0.00038/2 (exp.) +- 0.000120/2 (th.) = 0.00063

      call alfahat(mz,  dgamma,az)
      call alfahat(mw,  dgamma,aw)
      call alfahat(mb,  dgamma,ab)
      call alfahat(mtau,dgamma,atau)

         az =   az/pi1
         aw =   aw/pi1 
         ab =   ab/pi1
       atau = atau/pi1
         as = alphas(mc)/pi1
      astau = alphas(mtau)/pi1
        asb = alphas(mb)/pi1
        asw = alphas(mw)/pi1

      if (f4lqcd.eqv..false.) then
         as = as*(1 + cmatch*as**2)
      else
         as = as*(1 + cmatch*as**2 + c3m(3)*as**3)
      endif
      as = arun(3,mc,mtau,as*pi1)/pi1

      rmt2 = (mmu/mtau)**2

      etat = astau/(1 + 75*astau/76/atau)/4
      etab = asb  /(1 + 69*asb  /80/ab  )/4
      etaw = asw  /(1 + 69*asw  /17/aw  )/4

C  Leading order summation (Marciano and Sirlin, PRL 56, 22 1986):
C
C     sew = (1 - 38/ 9.d0*atau*dlog(mb/mtau))**(- 9/19.d0)*
C    .      (1 - 40/ 9.d0*ab  *dlog(mw/mb  ))**(- 9/20.d0)*
C    .      (1 - 17/36.d0*aw  *dlog(ratzw2 ))**(-36/17.d0)
C
C  Next-to-leading order summation (JE, 2002):

      sew = (1 - 38/ 9.d0*atau *dlog(mb/mtau))**(  9/19.d0*(etat - 1))*
     .      (1 - 40/ 9.d0*ab   *dlog(mw/mb  ))**(  9/20.d0*(etab - 1))*
     .      (1 - 17/36.d0*aw   *dlog(ratzw2 ))**( 36/17.d0*(etaw - 1))*
     .      (1 + 25/ 6.d0*astau*dlog(mb/mtau))**(- 9/19.d0*etat)*
     .      (1 + 23/ 6.d0*asb  *dlog(mw/mb)  )**(- 9/20.d0*etab)*
     .      (1 + 23/12.d0*asw  *dlog(ratzw2) )**(-36/17.d0*etaw)

C  Approximation:
C
C      sew = (1 - 38/ 9.d0*atau *dlog(mb/mtau))**(- 9/19.d0)*
C     .      (1 - 40/ 9.d0*ab   *dlog(mw/mb  ))**(- 9/20.d0)*
C     .      (1 - 17/36.d0*aw   *dlog(ratzw2 ))**(-36/17.d0)*
C     .      (1 + 25/ 6.d0*astau*dlog(mb/mtau))**(- 3/25.d0*atau)*
C     .      (1 + 23/ 6.d0*asb  *dlog(mz/mb)  )**(- 3/23.d0*ab)

      call cipqcd(as,delqcd,pterr)

       g0 = gf**2*mtau**5/192/pi1**3*(1 + 3*(mtau/mw)**2/5)

C Kai016 ---------------------------------------
C     Corrections to the W mass
C ----------------------------------------------

      if( modtype.eq.1 ) then
      
        g0 = g0 + 
     - ((fitcph**2 - fits2b)*gf**2*kkcc*mtau**5*
     -    ((gf*kkss*mtau**2)/atau)**0.2)/
     -  (160.*2**0.9*fitx*(-kkcc + kkss)*pi1**3.2)

      endif
      
      if( modtype.eq.2 ) then
      
        g0 = g0 + 
     -  ((fitcph**2 - 2*fits2b)*gf**2*kkcc*mtau**5*
     -    ((gf*kkss*mtau**2)/atau)**0.2)/
     -  (640.*2**0.9*fitx*(-kkcc + kkss)*pi1**3.2)

      endif

      if( modtype.eq.3 ) then

        g0 = g0 + 
     - (gf**2*((-1 + fitcph)**2 + 
     -      (fitsph**2*kkcc)/(-kkcc + kkss))*mtau**5*
     -    ((gf*kkss*mtau**2)/atau)**0.2)/
     -  (160.*2**0.9*fitx*pi1**3.2)

      endif

      if( modtype.eq.4 ) then

        g0 = g0 + 
     - (gf**2*((-1 + fitcph)**2 + 
     -      (fitsph**2*kkcc)/(-kkcc + kkss))*mtau**5*
     -    ((gf*kkss*mtau**2)/atau)**0.2)/
     -  (160.*2**0.9*fitx*pi1**3.2)

      endif

C Kai016 ---------------------------------------

       ge = g0*(1 + atau/2*(25/4.d0 - pi2))*(1 - 8*(me/mtau)**2)
      gmu = g0*(1 - 8*rmt2 + 8*rmt2**3 - rmt2**4 - 12*rmt2**2*dlog(rmt2)
     .    + atau/2*(25/4.d0 - pi2 - (68 + 24*dlog(rmt2))*rmt2
     .    + 32*pi2*(sqrt(rmt2)**3 + rmt2) 
     .    - (16*pi2 + 273 - 36*dlog(rmt2)*(1 - dlog(rmt2)))*rmt2**2))
      gud = 3*vud**2*g0*sew*(1 + atau/2*(85/12.d0 - pi2) + delqcd
c    .    - 0.00044d0  
c    .    - 0.10d0*(mu**2 + md**2)/(ms**2 - md**2)
     .    + 11*pi2/4*cdim4*(as/mtau**2)**2 - 24*pi2*cdim6 - 16*pi2*cdim8
     .    - 16*pi2*(fpipm*mpipm/mtau**2)**2*(1 + 23/8.d0*as**2)
     .    +  8*pi2*(fkpm*mkpm*as/mtau**2)**2)
c    .    - 0.0048d0)

C  The Delta S = - 1 contribution, rustot = 0.0285 +- 0.0007, is taken from 
C  experiment to avoid use of a poorly convergent series proportional to ms^2.

      gtot = (ge + gmu + gud)/(1 - rustot)

      taulif = hbar/gtot*1.d15
      therr2 = (hbar*3*vud**2*g0*sew*pterr/gtot**2*1.d15)**2 + remerr**2

      return
      end


      subroutine cipqcd(as,delqcd,pterr)

      implicit none
      integer i,j
      dimension aan(4)
      double precision as,rtc2,rtc4,rtc6,rst2,lntc2,lntc22,delqcd,pterr
      double precision cum0a2,cum0a3,cum0a4,cum0a5,cumca2,cumca3,cumca4
      double precision cumsa2,aan,msrun,mcrun,mbrun,rtb2,rtb4,rtb6,lntb2
      double precision cipter,fopter
      include 'common.f'

      rst2 = (msrun(mtau)/mtau)**2

      rtc2 = (mtau/mcrun(mtau))**2
      rtc4 = rtc2**2
      rtc6 = rtc2**3

      rtb2 = (mtau/mbrun(mtau))**2
      rtb4 = rtb2**2
      rtb6 = rtb2**3

      lntc2  = dlog(rtc2)
      lntb2  = dlog(rtb2)
      lntc22 = lntc2**2

      cum0a2 =   299/24.d0  -   9*zeta3
      cum0a3 = 58057/288.d0 - 779*zeta3/4 + 75*zeta5/2
      cum0a4 = 78631453/20736.d0 - 1704247*zeta3/432.d0 + 
     .     4185*zeta3**2/8.d0 + 34165*zeta5/96.d0 - 1995*zeta7/16.d0

      cum0a5 = 77.d0*3 ! needs better estimate

      cumca2 =       (107/4500.d0 - lntc2/225   )*rtc2
     .       -   (1597/7938000.d0 - lntc2/18900 )*rtc4
     .       + (3991/750141000.d0 - lntc2/595350)*rtc6
     .       +       (107/4500.d0 - lntb2/225   )*rtb2
     .       -   (1597/7938000.d0 - lntb2/18900 )*rtb4
     .       + (3991/750141000.d0 - lntb2/595350)*rtb6

      cumca3 = 
     . - (5.71561d-2 - 232933*lntc2/2916000   +  41*lntc22/3240   )*rtc2
     . - (9.96676d-4 +  57811*lntc2/342921600 - 187*lntc22/2721600)*rtc4
     . + (3.65225d-5 - 8.94071d-5*lntc2     + 361*lntc22/214326000)*rtc6

      cumsa2 = (32 - 24*zeta3)*rst2

C  exact coefficients (for 4 quark flavors) extracted from 
C  S. A. Larin, T. van Ritbergen and J. A. M. Vermaseren, NIKHEF-H/94-30:
C
C -        978751/34992000        +     zeta2/324       -   1847*zeta3/64800 
C -     112863433/17146080000     -  23*zeta2/272160    + 103693*zeta3/21772800
C   1097479447501/680527915200000 + 179*zeta2/107163000 -   4943*zeta3/3763200
C - 19315573*lntc2/2160406080000

      j = 1.d4*as
      if (j.lt.0) then
         do 20 i = 1, 4
            aan(i) = 0.d0
 20      continue
      else if (j.ge.2000) then
         do 30 i = 1, 4
            aan(i) = an(i,2000)
 30      continue
      else
         do 40 i = 1, 4
            aan(i) = an(i,j)+(as - asgrid(j))*(an(i,j+1) - an(i,j))*1.d4
 40      continue
      endif

c      delqcd =  aan(1) + cum0a2 *aan(2) + cum0a3*aan(3) + cum0a4*aan(4)

c      cipter = cum0a4*aan(4) ! needs calculation of aan(5)
c      pterr = cipter

C Fixed-order perturbation theory (FOPT):

      delqcd = as + 5.202d0*as**2 + 26.366d0*as**3 + 127.079d0*as**4     

      fopter = 127.079d0*as**4
      pterr = fopter

C "Pade-resummed" FOPT:

c      delqcd = as/(1 - 5.1*as) + 0.102*as**2 + 0.356*as**3 - 5.572*as**4 

c      delqcd = as/(1 - 5.202*as) - 0.695*as**3 - 13.691*as**4 

c      fopter = 127.079d0*as**4
c      pterr = fopter

C Weighted average of FOPT-CIPT:

c      pterr = 1.d0/dsqrt(1.d0/fopter**2 + 1.d0/cipter**2)

c      delqcd = ((as + 5.202d0*as**2 + 26.366d0*as**3 + 127.079d0*as**4)
c     .     /fopter**2
c     .     + (aan(1) + cum0a2 *aan(2) + cum0a3*aan(3) + cum0a4*aan(4))
c     .     /cipter**2)*pterr**2

      delqcd = delqcd + (cumca2 + cumsa2)*as**2  + cumca3*as**3

      return
      end
