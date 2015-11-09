      subroutine nue(q2,gvnue,ganue)

C   computes nu_mu electron scattering (at present only for Q^2 = 0).
C
C   following W. J. Marciano and A. Sirlin, PRD 22 (1980) 2695.
C         and W. J. Marciano in "Precision Tests of the SM".

      implicit none
      double precision q2,rhonue,s2tnue,sin4th,al,alh,cz,cgamma
      double precision gvnue,ganue,epse_L,epse_R,shat
      include 'common.f'

      if (q2.ne.0.d0) then
         print*,'presently nu-e scattering rad. corr. only at Q^2 = 0'
         return
      endif

          al = alpha/pi1
         alh = alphat/pi1
      sin4th = sinhat*sinhat

C   box contributions:      

          cz = 19/8.d0 -  7*sinhat/2 + 3*sin4th
      cgamma = 19/8.d0 - 17*sinhat/4 + 3*sin4th
      
      rhonue =  rhonc + alh/2/sinhat/coshat*cz
      s2tnue = shat(0.d0) - alh/2/coshat*cgamma

C   muon-neutrino charge-radius contribution:

      s2tnue = s2tnue - al/6*(dlog(mw2/mmu**2) + 1/6.d0) - 2*al/9

C   effective 4-Fermi neutral current interactions:

      gvnue = rhonue*( - 1/2.d0 + 2*s2tnue)
      ganue = rhonue*( - 1/2.d0)

C Kai018 ---------------------------------------
C     Corrections to the e couplings
C ----------------------------------------------

      if( modtype.eq.1 ) then

       gvnue = gvnue +
     - (kkcc*(i3(1) - fits2b*i3(1) + i3r(1) + 
     -       2*(-1 + fitcph)*(1 + fitcph*kkss)*q(1)) + 
     -    kkss*((-1 + fits2b)*i3(1) - i3r(1) + 
     -       2*(1 + fitcph*(-1 + kkss) - fits2b*kkss)*q(1)))/
     -  (fitx*(kkcc - kkss))

       ganue = ganue +
     - (i3(1) - fits2b*i3(1) - i3r(1))/fitx

      endif

      if( modtype.eq.2 ) then

       gvnue = gvnue +
     - (kkcc*(i3(1) - 2*fits2b*i3(1) + i3r(1) + 
     -       2*(-1 + fitcph)*(1 + fitcph*kkss)*q(1)) + 
     -    kkss*(-i3(1) + 2*fits2b*i3(1) - i3r(1) + 
     -       2*(1 + fitcph*(-1 + kkss) - 2*fits2b*kkss)*q(1)))/
     -  (4.*fitx*(kkcc - kkss))

       ganue = ganue +
     - (i3(1) - 2*fits2b*i3(1) - i3r(1))/(4.*fitx)

      endif

      if( modtype.eq.3 ) then
        
       gvnue = gvnue +
     - (-2*(-1 + fitcph)**2*kkcc*kkss)/(fitx*(kkcc - kkss))

       ganue = ganue + 0.d0

      endif

      if( modtype.eq.4 ) then

       gvnue = gvnue + ((kkcc - kkss)*i3(1) + 
     -    2*(-1 + fitcph)*(fitcph*kkcc - kkss)*kkss*q(1))/
     -  (fitx*(kkcc - kkss))

       ganue = ganue + i3(1)/fitx

      endif

C Kai018 ---------------------------------------

C   for PDG Table 10.3:
C
C      write(7,*),'rho_nue        = ',rhonue
C      write(7,*),'kappa_nue      = ',s2tnue/sinhat

      if (flagzp.ne.0) then
         epse_L = (gvnue + ganue)/2
         epse_R = (gvnue - ganue)/2
         epse_L = rhoeff*epse_L + 2*rhoezp*eps2_L(0)*eps2_L(1)
     .          + rhozzp*(eps2_L(1) + 2*eps2_L(0)*epse_L)
         epse_R = rhoeff*epse_R + 2*rhoezp*eps2_L(0)*eps2_R(1)
     .          + rhozzp*(eps2_R(1) + 2*eps2_L(0)*epse_R)
          gvnue = epse_L + epse_R
          ganue = epse_L - epse_R
      endif

      return
      end

