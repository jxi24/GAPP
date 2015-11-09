      subroutine apv(Qw,Z,AA,C1u,C1d,C2u,C2d)

C   following W. J. Marciano and A. Sirlin, PRD 29 (1984) 75
C         and W. J. Marciano, TASI 1993 lectures.

      implicit none
      integer Z,AA
      double precision Qw,C1u,C1d,C2u,C2d,rhoapv,s2tapv,alphas,asw,asz
      double precision al,alh,cgzbx1,cgzbx2,gzboxh,lamd2u,lamd2d,m2,qwp
      double precision qwn,shat
      include 'common.f'

      al  =  alpha/pi1
      alh = alphat/pi1
      asw = alphas(mw)/pi1
      asz = alphas(mz)/pi1

      cgzbx1 = 0.29418d0
      cgzbx2 =     11.d0
C     cgzbx2 =      8.693d0   ! for PVDIS-JLab
      gzboxh = alphat/pi1*(1 - 4*sinhat)*cgzbx1

C   for PDG Table 10.3:
C
C      rhoapv = rhonc - al/2 - alh/2*(1/sinhat + 4*(1 - 4*sinhat)*cgzbx2
C     .       + 9/coshat/8*(1 - 16*sinhat/9)*(1/sinhat - 4 + 8*sinhat))
C
C      s2tapv = shat(0.d0) - 2*al/9 + alh/4*(2 - 9/sinhat/4
C     .       + (1 - 4*sinhat)/3*(dlog(mz2/me**2) + 1/6.d0)
C     .       - (9/2.d0 - 8*sinhat)*(1 - 4*sinhat)*cgzbx2
C     .       - 9/coshat/8*(1/sinhat - 4 + 8*sinhat)
C     .       * (1 - 4*sinhat + 32*sinhat**2/9))

C   for C2u and C2d:

          m2 = 87.d-3**2
C         m2 = 6.5d0         ! for PVDIS-JLab
      lamd2u = - al/9 *(1 - 4*sinhat)
      lamd2d =   al/36*(1 - 4*sinhat)

C   box contributions:

      lamd2u = lamd2u + al/2*(1/sinhat + (1 - 8*sinhat/3)*cgzbx2 
     .       + (1 - 4*sinhat)/coshat/16*(3/sinhat - 8 + 32*sinhat/3))
      lamd2d = lamd2d + al/2*( - 1/sinhat/4 + (1 - 4*sinhat/3)/2*cgzbx2 
     .       + (1 - 4*sinhat)/coshat/16*(3/sinhat - 4 + 8*sinhat/3))

C   charge-radius contributions:

      lamd2u = lamd2u +    al/9*((dlog(mw2/m2) + 1/6.d0) 
     .       - (1 - 8*sinhat/3)*( dlog(mz2/m2) + 1/6.d0))
      lamd2d = lamd2d -  2*al/9*((dlog(mw2/m2) + 1/6.d0)
     .       - (1 - 4*sinhat/3)*( dlog(mz2/m2) + 1/6.d0)/4)

C   for PDG Table 10.3:
C      
C      write(7,*),'rho^prime      = ',rhoapv
C      write(7,*),'rho            = ',rhonc
C      write(7,*),'kappa^prime    = ',s2tapv/sinhat
C      write(7,*),'kappa          = ',(shat(0.d0) - 2*al/9)/sinhat
C      write(7,*),'lambda_1d      = ',2*gzboxh/3
C      write(7,*),'lambda_2u      = ', - lamd2u
C      write(7,*),'lambda_2d      = ', - lamd2d
C      write(7,*)

C     C1u = -    (1/2.d0*rhoapv*(1 - 8/3.d0*s2tapv) +   gzboxh/3)
C     C1d = - ( - 1/2.d0*rhoapv*(1 - 4/3.d0*s2tapv) - 2*gzboxh/3)
      C2u = -    (1/2.d0*rhonc *(1 - 4*shat(0.d0) + 8*al/9) + lamd2u)
      C2d = - ( - 1/2.d0*rhonc *(1 - 4*shat(0.d0) + 8*al/9) + lamd2d)
      C2u = -    (1/2.d0*rhonc *(1 - 4*shat(0.d0) + 8*al/9) + lamd2u)
      C2d = - ( - 1/2.d0*rhonc *(1 - 4*shat(0.d0) + 8*al/9) + lamd2d)

      qwp = (rhonc - al/2)*(1 - 4*shat(0.d0) + 8*al/9
     .    - al/3*(1 - 4*sinhat)*(dlog(mz2/me**2) + 1/6.d0))
     .    + alh/4/sinhat*((2 + 5*(1 - asw)) + (9/4.d0 - 5*sinhat)*(1 
     .    - 4*sinhat + 8*sinhat**2)/coshat*(1 - asz))
     .    + 5*alh/2*(1 - 4*sinhat)*cgzbx2

      qwn = - (rhonc - al/2 - alh/2*(1/sinhat*(1 - 2*asw) 
     .      + 4*(1 - 4*sinhat)*cgzbx2 + 9/coshat/8*(1 - 16*sinhat/9)*
     .        (1/sinhat - 4 + 8*sinhat)*(1 - asz)) + gzboxh)

      C1u = - (qwp - qwn/2)/3
      C1d = - (qwn - qwp/2)/3

C Kai019 ---------------------------------------
C     Corrections to the e couplings
C ----------------------------------------------

       if( modtype.eq.1 ) then

       C1u = C1u +         
     - (2*(kkcc*(-(i3r(1)*(i3(4) + i3r(4))) + 
     -         i3(1)*(i3(4) - fits2b*i3(4) + i3r(4)) + 
     -         2*((-1 + fitcph)*(1 + fitcph*kkss)*i3(1) + 
     -            i3r(1) + fitcph*(-1 + kkss)*i3r(1))*q(4)) + 
     -      kkss*((-1 + fits2b)*i3(1)*i3(4) - i3(1)*i3r(4) + 
     -         i3r(1)*(i3(4) + i3r(4)) + 
     -         2*((1 + fitcph*(-1 + kkss) - fits2b*kkss)*
     -             i3(1) + (-1 + fitcph - fitcph*kkss)*i3r(1))*
     -          q(4))))/(fitx*(kkcc - kkss))

       C1d = C1d +
     - (2*(kkcc*(-(i3r(1)*(i3(7) + i3r(7))) + 
     -         i3(1)*(i3(7) - fits2b*i3(7) + i3r(7)) + 
     -         2*((-1 + fitcph)*(1 + fitcph*kkss)*i3(1) + 
     -            i3r(1) + fitcph*(-1 + kkss)*i3r(1))*q(7)) + 
     -      kkss*((-1 + fits2b)*i3(1)*i3(7) - i3(1)*i3r(7) + 
     -         i3r(1)*(i3(7) + i3r(7)) + 
     -         2*((1 + fitcph*(-1 + kkss) - fits2b*kkss)*
     -             i3(1) + (-1 + fitcph - fitcph*kkss)*i3r(1))*
     -          q(7))))/(fitx*(kkcc - kkss))

       C2u = C2u +
     - (2*(kkcc - kkss)*i3r(4)*
     -     (-i3(1) - i3r(1) + 2*(1 + fitcph*(-1 + kkss))*q(1))
     -     + 2*i3(4)*(kkcc*
     -        (i3(1) - fits2b*i3(1) + i3r(1) + 
     -          2*(-1 + fitcph)*(1 + fitcph*kkss)*q(1)) + 
     -       kkss*((-1 + fits2b)*i3(1) - i3r(1) + 
     -          2*(1 + fitcph*(-1 + kkss) - fits2b*kkss)*q(1)))
     -    )/(fitx*(kkcc - kkss))
     
       C2d = C2d +
     - (2*(kkcc - kkss)*i3r(7)*
     -     (-i3(1) - i3r(1) + 2*(1 + fitcph*(-1 + kkss))*q(1))
     -     + 2*i3(7)*(kkcc*
     -        (i3(1) - fits2b*i3(1) + i3r(1) + 
     -          2*(-1 + fitcph)*(1 + fitcph*kkss)*q(1)) + 
     -       kkss*((-1 + fits2b)*i3(1) - i3r(1) + 
     -          2*(1 + fitcph*(-1 + kkss) - fits2b*kkss)*q(1)))
     -    )/(fitx*(kkcc - kkss))

      endif

      if( modtype.eq.2 ) then

       C1u = C1u +
     - (kkcc*(-(i3r(1)*(i3(4) + i3r(4))) + 
     -       i3(1)*(i3(4) - 2*fits2b*i3(4) + i3r(4)) + 
     -       2*((-1 + fitcph)*(1 + fitcph*kkss)*i3(1) + 
     -          i3r(1) + fitcph*(-1 + kkss)*i3r(1))*q(4)) + 
     -    kkss*(i3(1)*(-i3(4) + 2*fits2b*i3(4) - i3r(4)) + 
     -       i3r(1)*(i3(4) + i3r(4)) + 
     -       2*((1 + fitcph*(-1 + kkss) - 2*fits2b*kkss)*
     -           i3(1) + (-1 + fitcph - fitcph*kkss)*i3r(1))*
     -        q(4)))/(2.*fitx*(kkcc - kkss))

       C1d = C1d +  
     - (kkcc*(-(i3r(1)*(i3(7) + i3r(7))) + 
     -       i3(1)*(i3(7) - 2*fits2b*i3(7) + i3r(7)) + 
     -       2*((-1 + fitcph)*(1 + fitcph*kkss)*i3(1) + 
     -          i3r(1) + fitcph*(-1 + kkss)*i3r(1))*q(7)) + 
     -    kkss*(i3(1)*(-i3(7) + 2*fits2b*i3(7) - i3r(7)) + 
     -       i3r(1)*(i3(7) + i3r(7)) + 
     -       2*((1 + fitcph*(-1 + kkss) - 2*fits2b*kkss)*
     -           i3(1) + (-1 + fitcph - fitcph*kkss)*i3r(1))*
     -        q(7)))/(2.*fitx*(kkcc - kkss))

       C2u = C2u +  
     - ((kkcc - kkss)*i3r(4)*
     -     (-i3(1) - i3r(1) + 2*(1 + fitcph*(-1 + kkss))*q(1))
     -     + i3(4)*(kkcc*(i3(1) - 2*fits2b*i3(1) + i3r(1) + 
     -          2*(-1 + fitcph)*(1 + fitcph*kkss)*q(1)) + 
     -       kkss*(-i3(1) + 2*fits2b*i3(1) - i3r(1) + 
     -          2*(1 + fitcph*(-1 + kkss) - 2*fits2b*kkss)*q(1)
     -          )))/(2.*fitx*(kkcc - kkss))

       C2d = C2d +  
     - ((kkcc - kkss)*i3r(7)*
     -     (-i3(1) - i3r(1) + 2*(1 + fitcph*(-1 + kkss))*q(1))
     -     + i3(7)*(kkcc*(i3(1) - 2*fits2b*i3(1) + i3r(1) + 
     -          2*(-1 + fitcph)*(1 + fitcph*kkss)*q(1)) + 
     -       kkss*(-i3(1) + 2*fits2b*i3(1) - i3r(1) + 
     -          2*(1 + fitcph*(-1 + kkss) - 2*fits2b*kkss)*q(1)
     -          )))/(2.*fitx*(kkcc - kkss))

      endif

      if( modtype.eq.3 ) then
       
       C1u = C1u +         
     - (4*(-1 + fitcph)**2*kkcc*kkss*i3(1)*q(4))/
     -  (fitx*(kkcc - kkss))

       C1d = C1d +
     - (4*(-1 + fitcph)**2*kkcc*kkss*i3(1)*q(7))/
     -  (fitx*(kkcc - kkss))

       C2u = C2u +
     - (4*(-1 + fitcph)*kkss*(-(fitcph*kkcc) + kkss)*i3(4)*
     -    q(1))/(fitx*(-kkcc + kkss))
     
       C2d = C2d +
     - (4*(-1 + fitcph)*kkss*(-(fitcph*kkcc) + kkss)*i3(7)*
     -    q(1))/(fitx*(-kkcc + kkss))

      endif

      if( modtype.eq.4 ) then
       
       C1u = C1u +         
     - (2*(kkcc - kkss)*i3(1)*i3(4) + 
     -    4*(-1 + fitcph)*(fitcph*kkcc - kkss)*kkss*i3(1)*q(4)
     -    )/(fitx*(kkcc - kkss))

       C1d = C1d +
     - (2*(kkcc - kkss)*i3(1)*i3(7) + 
     -    4*(-1 + fitcph)*(fitcph*kkcc - kkss)*kkss*i3(1)*q(7)
     -    )/(fitx*(kkcc - kkss))

       C2u = C2u +
     - (2*i3(4)*((kkcc - kkss)*i3(1) + 
     -      2*(-1 + fitcph)*(fitcph*kkcc - kkss)*kkss*q(1)))/
     -  (fitx*(kkcc - kkss))
     
       C2d = C2d +
     - (2*i3(7)*((kkcc - kkss)*i3(1) + 
     -      2*(-1 + fitcph)*(fitcph*kkcc - kkss)*kkss*q(1)))/
     -  (fitx*(kkcc - kkss))

      endif

C Kai019 ---------------------------------------

      if (flagzp.ne.0) then
         C1u = rhoeff*C1u + 2*rhoezp*a2(1)*v2(4)
     .       - rhozzp*(v2(4) - a2(1)*(1 - 8/3.d0*shat(0.d0)))
         C1d = rhoeff*C1d + 2*rhoezp*a2(1)*v2(7)
     .       - rhozzp*(v2(7) + a2(1)*(1 - 4/3.d0*shat(0.d0)))
         C2u = rhoeff*C2u + 2*rhoezp*a2(4)*v2(1)
     .       + rhozzp*(v2(1) - a2(4)*(1 - 4*shat(0.d0)))
         C2d = rhoeff*C2d + 2*rhoezp*a2(7)*v2(1)
     .       - rhozzp*(v2(1) + a2(7)*(1 - 4*shat(0.d0)))
      endif
 
      Qw = - 2*((AA + Z)*C1u + (2*AA - Z)*C1d)

      return
      end


      subroutine moller(qweake,q2,y)

C   rad. corr. following A. Czarnecki and W. J. Marciano, preprint TTP95-14.

      implicit none
      double precision alh,f1half,s2tmol,alrmol,q2,qweake
      double precision y,y2,y3,y4,y5,log1,log2,fy,shat,alfaq2,dgamma
      include 'common.f'
      
C     SLAC-E-158 parameters:
C
C     q2 = y*2*me*46.5d0 (approximately)
C     q2 = 0.026d0 (MC) 
C     y  = 0.599d0 (adjusted to reproduce MC analyzing power = 3.25 ppm)
C
C     JLab Moller parameters:      
C
C     q2 = y*2*me*11.d0 (approximately)
C     q2 = 0.0056d0 (MC)
C     y  = 0.571d0 (adjusted to reproduce MC asymmetry = 35.6 ppb
C                   corresponding to analyzing power of 759 ppb)

      call alfahat(dsqrt(q2),dgamma,alfaq2)
      alh = alphat/pi1
      
      y2 = y*y
      y3 = y*y2
      y4 = y*y3
      y5 = y*y4
      log1 = dlog(y)
      log2 = dlog(1 - y)

      f1half = 17*zeta2/2 + 70*dlog(2.d0)/9 - 8*dlog(2.d0)**2/3 ! (max. asym.)
      
      fy = - 2/3.d0*(log1 + log2) + (
     .     - 2*(1 - y)*(3 - 3*y         +  4*y3 -  3*y4       )*log2
     .     -       2*y*(1 + 3*y -  6*y2 +  8*y3 -  3*y4       )*log1
     .     +   (1 - y)*(2 - 2*y -  7*y2 + 10*y3 -  8*y4 + 3*y5)*log2**2
     .     -         y*(2 - 3*y -  5*y2 +  8*y3 -  7*y4 + 3*y5)*log1**2
     .     +           (2 - 4*y         + 11*y3 - 13*y4 + 9*y5 - 3*y*y5)
     .                 *(pi2 - 2*log1*log2))/(1 - y + y2)**2

      s2tmol = shat(0.d0) - 2*alpha/9/pi1 - alh/16*(1/sinhat
     .       - (1 - 4*shat(0.d0))*(3/coshat/4*(1/sinhat - 4 + 8*sinhat)
     .       + (22*dlog(y*mz2/q2)/3 + 85/9.d0 + fy)*alpha/alphat))

C Alternative (to solve one needs resummation or 2-loop result):
C
C      s2tmol = shat(0.d0) - 2*alpha/9/pi1 - alh/16*(1/sinhat
C     .       - (1 - 4*sinhat)*(3/coshat/4*(1/sinhat - 4 + 8*sinhat)
C     .       + (22*dlog(y*mz2/q2)/3 + 85/9.d0 + fy)))
C
C Alternative (suggestion by Czarnecki and Marciano):
C
C      s2tmol = shat(0.d0) - 2*alpha/9/pi1 - alh/16*(1/sinhat
C     .       - (1 - 2*shat(0.d0) - 2*sinhat)
C     .     *(3/coshat/4*(1/sinhat - 4 + 8*sinhat)
C     .       + (22*dlog(y*mz2/q2)/3 + 85/9.d0 + fy)*alpha/alphat))

      qweake = - rhonc*(1 - 4*s2tmol)

C Kai025 ---------------------------------------
C     Corrections to the weak charge of the electron
C----------------------------------------------

      if( modtype.eq.1 ) then

         qweake = qweake + 
     - (4*(kkcc - kkss)*
     -     ((-1 + fits2b)*i3(1)**2 + i3r(1)**2) - 
     -    8*(((-1 + fitcph)*kkcc*(1 + fitcph*kkss) + 
     -          kkss*(1 + fitcph*(-1 + kkss) - fits2b*kkss))*
     -        i3(1) + (1 + fitcph*(-1 + kkss))*(kkcc - kkss)*
     -        i3r(1))*q(1))/(fitx*(kkcc - kkss))

      endif

      if( modtype.eq.2 ) then

         qweake = qweake +
     -  ((kkcc - kkss)*((-1 + 2*fits2b)*i3(1)**2 + 
     -       i3r(1)**2) - 
     -    2*(((-1 + fitcph)*kkcc*(1 + fitcph*kkss) + 
     -          kkss*(1 + fitcph*(-1 + kkss) - 2*fits2b*kkss))*
     -        i3(1) + (1 + fitcph*(-1 + kkss))*(kkcc - kkss)*
     -        i3r(1))*q(1))/(fitx*(kkcc - kkss))

      endif

      if( modtype.eq.3 ) then

         qweake = qweake +
     -  (-8*(-1 + fitcph)**2*kkcc*kkss*i3(1)*q(1))/
     -  (fitx*(kkcc - kkss))

      endif

      if( modtype.eq.4 ) then

         qweake = qweake +
     -  (4*i3(1)*((-kkcc + kkss)*i3(1) + 
     -      2*(-1 + fitcph)*kkss*(-(fitcph*kkcc) + kkss)*q(1))
     -    )/(fitx*(kkcc - kkss))

      endif
        
C ----------------------------------------------


      alrmol = qweake*gf*q2/dsqrt(2.d0)/pi1/alfaq2
     .       * (1 - y)/(1 + y4 + (1 - y)**4)*1.01d0

      if (flagzp.ne.0) qweake = rhoeff*qweake - 4*rhoezp*a2(1)*v2(1)
     .     + 2*rhozzp*(v2(1) + a2(1)*(1 - 4.d0*shat(0.d0)))

      return
      end
