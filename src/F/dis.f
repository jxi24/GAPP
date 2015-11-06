      subroutine nuhnutev(NuTeV,epsu_L,epsd_L,epsu_R,epsd_R)

C   momentum transfer and weight factors for NuTeV experiment as in
C   hep-ex/9806013:

      implicit none
      double precision NuTeV,q2NuTV,cuL2,cuR2,cdL2,cdR2
      double precision epsu_L,epsd_L,epsu_R,epsd_R

      q2NuTV = - 20.d0
        cuL2 =   0.8587d0
        cdL2 =   0.8828d0
        cuR2 = - 1.1657d0
        cdR2 = - 1.2288d0

      call nuh(q2NuTV,epsu_L,epsd_L,epsu_R,epsd_R)

      NuTeV = cuL2*epsu_L**2 + cdL2*epsd_L**2 
     .      + cuR2*epsu_R**2 + cdR2*epsd_R**2

      return
      end


      subroutine nuhccfr(CCFR)

C   momentum transfer and weight factors for CCFR experiment as in
C   FNAL-Pub-97/001-E:

      implicit none
      double precision CCFR,q2CCFR,cgL2,cgR2,cdelL2,cdelR2
      double precision epsu_L,epsd_L,epsu_R,epsd_R,gL2,gR2,delL2,delR2

      q2CCFR = - 35.d0
        cgL2 =   1.7897d0
        cgR2 =   1.1479d0
      cdelL2 = - 0.0916d0
      cdelR2 = - 0.0782d0

      call nuh(q2CCFR,epsu_L,epsd_L,epsu_R,epsd_R)

         gL2 = epsu_L**2 + epsd_L**2
         gR2 = epsu_R**2 + epsd_R**2
       delL2 = epsu_L**2 - epsd_L**2
       delR2 = epsu_R**2 - epsd_R**2

        CCFR = cgL2*gL2 + cgR2*gR2 + cdelL2*delL2 + cdelR2*delR2

      return
      end


      subroutine nuhcdhs(Rnu,Rnubar,Rnuba2)

C   average momentum transfers as given by F. Perrier in 
C   "Precision Tests of the Standard Electroweak Model",
C   ed. P. Langacker (World Scientific, 1995) p. 466 
C
C   weight factors for CDHS experiment as in ZPC 45 (1990) 361.

      implicit none
      double precision q2,Rnu,Rnubar,Rnuba2,delta,a_Lu,a_Ld,a_Ru,a_Rd
      double precision epsu_L,epsd_L,epsu_R,epsd_R

         q2 = - 21.d0
      delta = 0.023d0
       a_Lu = 0.936d0
       a_Ld = 1.045d0
       a_Ru = 0.379d0
       a_Rd = 0.453d0

      call nuh(q2,epsu_L,epsd_L,epsu_R,epsd_R)

      Rnu = (1.d0 - delta)*(a_Lu*epsu_L**2 + a_Ld*epsd_L**2 
     .                    + a_Ru*epsu_R**2 + a_Rd*epsd_R**2)

         q2 = - 11.d0
      delta = 0.026d0
       a_Lu = 0.948d0
       a_Ld = 1.134d0
       a_Ru = 2.411d0
       a_Rd = 2.690d0

      call nuh(q2,epsu_L,epsd_L,epsu_R,epsd_R)

      Rnubar = (1.d0 - delta)*(a_Lu*epsu_L**2 + a_Ld*epsd_L**2 
     .                       + a_Ru*epsu_R**2 + a_Rd*epsd_R**2)

         q2 = - 11.d0
      delta = 0.024d0
       a_Lu = 0.944d0
       a_Ld = 1.126d0
       a_Ru = 2.295d0
       a_Rd = 2.563d0

      call nuh(q2,epsu_L,epsd_L,epsu_R,epsd_R)

      Rnuba2 = (1.d0 - delta)*(a_Lu*epsu_L**2 + a_Ld*epsd_L**2 
     .                       + a_Ru*epsu_R**2 + a_Rd*epsd_R**2)

      return
      end


      subroutine nuh(q2,epsu_L,epsd_L,epsu_R,epsd_R)

C   following W. J. Marciano and A. Sirlin, PRD 22 (1980) 2695.

      implicit none
      double precision q2,rhonuh,s2efq2,sineff,sin4th,sin6th
      double precision dgamma,alfaq2,alq,alh,az,agamma,abetal,abetar
      double precision epsu_L,epsd_L,epsu_R,epsd_R
      include 'common.f'

      alfaq2 = alphat
      call alfahat(dsqrt(-q2),dgamma,alfaq2)

         alh = alphat/pi1
         alq = alfaq2/pi1
      sin4th = sinhat*sinhat
      sin6th = sinhat*sin4th

C   box contributions:      

          az = (5/2.d0 - 15*sinhat/4 - sin4th/5 + 14*sin6th/9)/2/coshat
      agamma = (5/2.d0 -61*sinhat/20-9*sin4th/10+ 14*sin6th/9)/2/coshat
      abetal = - (9/16.d0 - 3*sinhat/4 + 4*sin4th/15)/coshat
      abetar = - 3*sin4th/10

      rhonuh = rhonc + alh/2/coshat/sinhat*az
      s2efq2 = sineff(q2,alq) - alh/2/coshat*agamma

C   electron-neutrino charge-radius contribution:

      s2efq2 = s2efq2 - alq/6*(dlog(-mw2/q2) + 8/3.d0)

C   effective 4-Fermi neutral current interactions:

      epsu_L = rhonuh*(   1/2.d0 - 2*s2efq2/3) 
     .       + alh/2/coshat/sinhat*abetal*(coshat + sinhat/3)
      epsd_L = rhonuh*( - 1/2.d0 +   s2efq2/3) 
     .       + alh/2/coshat/sinhat*abetal*(coshat - sinhat/3)
      epsu_R = - 2*rhonuh*s2efq2/3 - alh/6/coshat/sinhat*abetar
      epsd_R =     rhonuh*s2efq2/3 - alh/3/coshat/sinhat*abetar

C Kai017 ---------------------------------------
C     Corrections to the u and d couplings
C ----------------------------------------------

      if( modtype.eq.1 ) then

       epsu_L = epsu_L +
     - (-((-1 + fits2b)*(kkcc - kkss)*i3(4)) + 
     -    ((-1 + fitcph)*kkcc*(1 + fitcph*kkss) + 
     -       kkss*(1 + fitcph*(-1 + kkss) - fits2b*kkss))*q(4))
     -   /(fitx*(kkcc - kkss))

       epsu_R = epsu_R +
     - (kkcc*(i3r(4) + 
     -       (-1 + fitcph)*(1 + fitcph*kkss)*q(4)) + 
     -    kkss*(-i3r(4) + 
     -       (1 + fitcph*(-1 + kkss) - fits2b*kkss)*q(4)))/
     -  (fitx*(kkcc - kkss))

       epsd_L = epsd_L +
     - (-((-1 + fits2b)*(kkcc - kkss)*i3(7)) + 
     -    ((-1 + fitcph)*kkcc*(1 + fitcph*kkss) + 
     -       kkss*(1 + fitcph*(-1 + kkss) - fits2b*kkss))*q(7))
     -   /(fitx*(kkcc - kkss))

       epsd_R = epsd_R +
     - (kkcc*(i3r(7) + 
     -       (-1 + fitcph)*(1 + fitcph*kkss)*q(7)) + 
     -    kkss*(-i3r(7) + 
     -       (1 + fitcph*(-1 + kkss) - fits2b*kkss)*q(7)))/
     -  (fitx*(kkcc - kkss))
     
      endif
      
      if( modtype.eq.2 ) then
      
       epsu_L = epsu_L +
     - (-((-1 + 2*fits2b)*(kkcc - kkss)*i3(4)) + 
     -    ((-1 + fitcph)*kkcc*(1 + fitcph*kkss) + 
     -       kkss*(1 + fitcph*(-1 + kkss) - 2*fits2b*kkss))*
     -     q(4))/(4.*fitx*(kkcc - kkss))

       epsu_R = epsu_R +
     - (kkcc*(i3r(4) + 
     -       (-1 + fitcph)*(1 + fitcph*kkss)*q(4)) + 
     -    kkss*(-i3r(4) + 
     -       (1 + fitcph*(-1 + kkss) - 2*fits2b*kkss)*q(4)))/
     -  (4.*fitx*(kkcc - kkss))

       epsd_L = epsd_L +
     - (-((-1 + 2*fits2b)*(kkcc - kkss)*i3(7)) + 
     -    ((-1 + fitcph)*kkcc*(1 + fitcph*kkss) + 
     -       kkss*(1 + fitcph*(-1 + kkss) - 2*fits2b*kkss))*
     -     q(7))/(4.*fitx*(kkcc - kkss))

       epsd_R = epsd_R +
     - (kkcc*(i3r(7) + 
     -       (-1 + fitcph)*(1 + fitcph*kkss)*q(7)) + 
     -    kkss*(-i3r(7) + 
     -       (1 + fitcph*(-1 + kkss) - 2*fits2b*kkss)*q(7)))/
     -  (4.*fitx*(kkcc - kkss))
     
      endif


      if( modtype.eq.3 ) then

       epsu_L = epsu_L +
     - ((-1 + fitcph)**2*kkcc*kkss*q(4))/(fitx*(kkcc - kkss))

       epsu_R = epsu_R +
     - ((-1 + fitcph)**2*kkcc*kkss*q(4))/(fitx*(kkcc - kkss))

       epsd_L = epsd_L +
     - ((-1 + fitcph)**2*kkcc*kkss*q(7))/(fitx*(kkcc - kkss))

       epsd_R = epsd_R +
     - ((-1 + fitcph)**2*kkcc*kkss*q(7))/(fitx*(kkcc - kkss))
     
      endif

      if( modtype.eq.4 ) then

       epsu_L = epsu_L +
     - ((kkcc - kkss)*i3(4) + 
     -    (-1 + fitcph)*(fitcph*kkcc - kkss)*kkss*q(4))/
     -  (fitx*(kkcc - kkss))

       epsu_R = epsu_R +
     - ((-1 + fitcph)*(fitcph*kkcc - kkss)*kkss*q(4))/
     -  (fitx*(kkcc - kkss))

       epsd_L = epsd_L +
     - ((kkcc - kkss)*i3(7) + 
     -    (-1 + fitcph)*(fitcph*kkcc - kkss)*kkss*q(7))/
     -  (fitx*(kkcc - kkss))

       epsd_R = epsd_R +
     - ((-1 + fitcph)*(fitcph*kkcc - kkss)*kkss*q(7))/
     -  (fitx*(kkcc - kkss))

      endif

C Kai017 ---------------------------------------

      if (flagzp.ne.0) then
         epsu_L = rhoeff*epsu_L + 2*rhoezp*eps2_L(0)*eps2_L(4)
     .          + rhozzp*(eps2_L(4) + 2*eps2_L(0)*epsu_L)
         epsd_L = rhoeff*epsd_L + 2*rhoezp*eps2_L(0)*eps2_L(7)
     .          + rhozzp*(eps2_L(7) + 2*eps2_L(0)*epsd_L)
         epsu_R = rhoeff*epsu_R + 2*rhoezp*eps2_L(0)*eps2_R(4)
     .          + rhozzp*(eps2_R(4) + 2*eps2_L(0)*epsu_R)
         epsd_R = rhoeff*epsd_R + 2*rhoezp*eps2_L(0)*eps2_R(7)
     .          + rhozzp*(eps2_R(7) + 2*eps2_L(0)*epsd_R)
      endif

C   for PDG Table 10.3:
C
C      write(7,*),'rho_nuN        = ',rhonuh
C      write(7,*),'kappa_nuN      = ',s2efq2/(1-mw2/mz2)
C      write(7,*),'kappa_nuN (MS) = ',s2efq2/sinhat
C      write(7,*),'Q^2            = ',q2
C      write(7,*),'Lambda_uL      =
C     .     ',alh/2/coshat/sinhat*abetal*(coshat+sinhat/3)
C      write(7,*),'Lambda_dL      =
C     .     ',alh/2/coshat/sinhat*abetal*(coshat-sinhat/3)
C      write(7,*),'Lambda_uR = ', - alh/6/coshat/sinhat*abetar
C      write(7,*),'Lambda_dR = ', - alh/3/coshat/sinhat*abetar
C      write(7,*)

      return
      end


      double precision function sineff(q2,alq)

C  calculates effective weak angle in deep inelastic neutrino nucleon
C  scattering for ms^2 << -q2 << mt^2.

      implicit none
      integer f
      double precision q2,qq,alq,alh,alphas,aq,at,xa2mt0
      double precision mcrun,mbrun,rtauq2,ratcq2,ratbq2,ratwq2,rcs
      double precision const1,const2,lambda,pitv1,shat
      include 'common.f'

          qq = dsqrt(-q2)

         alh = alphat/pi1
          aq = alphas(qq)/pi1
          at = alphas(mt)/pi1
      const1 = (5/3.d0 - dlog(-q2/mz2))/3
      const2 = (55/12.d0 - 4*zeta3 - dlog(-q2/mz2))/3
      rtauq2 = - mtau**2/q2
      ratcq2 = - mcrun(qq)**2/q2
      ratbq2 = - mbrun(qq)**2/q2 
      ratwq2 =   mw2/q2
         rcs = dsqrt(coshat*sinhat)

      xa2mt0 = 1.d0
      if (fa2mt0.eqv..false.) xa2mt0 = 0.d0
      if (fasmt0.eqv..false.) then
         aq = 0.d0
         at = 0.d0
      endif

      sineff = sinhat
      do 100 f = 1, 9
         if (f.ne.6) then
            sineff = sineff + rcs*alq*nc(f)*q(f)*v(f)*const1
            if (f.ge.4) sineff = sineff + 3*rcs*alq*aq*q(f)*v(f)*const2
            if (flagmf.eqv..true.) then
               if (f.eq.3) sineff = sineff + rcs*alq*q(f)*v(f)*
     .                              (pitv1(rtauq2) - dlog(rtauq2))/3
               if (f.eq.5) sineff = sineff + rcs*alq*q(f)*v(f)*
     .                              (pitv1(ratcq2) - dlog(ratcq2))
               if (f.eq.9) sineff = sineff + rcs*alq*q(f)*v(f)*
     .                              (pitv1(ratbq2) - dlog(ratbq2))
            endif
         endif
 100  continue
      
      if ((flagmr.eqv..false.).or.(mt.le.mz))
     .     sineff = sineff - alq/6*(1 - 8*sinhat/3)*(dlog(rattz2)
     .     *(1 + at + alh*xa2mt0/3) - 13/12.d0*at - 5*alh*xa2mt0/4)

      sineff = sineff + alq/4*( - (3*coshat + 1/6.d0)*dlog(ratzw2)
     .                - 1/9.d0 + (lambda(1.d0/ratwq2) - 1)*
     .          (8*ratwq2*(coshat + 1/3.d0) + 6*coshat + 1/3.d0))

      sineff = sineff - alq*coshat*dlog(ratzw2)

      return
      end

      
      double precision function pitv1(x)

      implicit none
      double precision x,root

       root = dsqrt(1 + 4*x)
      pitv1 = - 4*x + root*(1 - 2*x)*dlog(4*x/(1 + root)**2)

      return
      end
