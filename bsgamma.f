      subroutine bsgamma(sgcenu)

C   computes R = B(b -> s gamma)/B(b -> c e nu)
C
C   following F.M. Borzumati and C. Greub, hep-ph/9802391
C          A. Czarnecki and W.J. Marciano, hep-ph/9804252
C               A.L. Kagan and M. Neubert, hep-ph/9805303

      implicit none
      integer i
      double precision sgcenu,alpem,ratckm,fphsp,hphsp,dabs2,dgamma
      double precision alphas,asb,eta,eta7(12),msrun,mcrun,mtrun
      double precision h3(8),h4(8),h5(8),h6(8),h7(8),h8(8),a7(8),e7(8)
      double precision f7(8),kl7(8),C70,C80,C41,C71,C81,C10eff,C20eff
      double precision C30eff,C40eff,C50eff,C60eff,C70eff,C80eff,C71eff
      double precision x,x2,x3,x4,x5,x12,x13,x14,x15,lnx,li2
      double precision z,z2,z3,z4,z5,z6,logz,logz2,logz3,l3p,rootz
      double precision f22,f221,f27,f271,f28,f77,f78,f88,cbue1,cbue2
      double precision Abrems,sum70,sum71,sum80,ftw,C7em,rtckm2
      double precision del,del2,del3,logdel,logdl2,li2del,sudfac
      double precision lamnp2,delnpg,delnpc,dlnpsl
      complex*16 spence,Vrtuel,r2,r7,r8
      include 'common.f'

      ratckm = 0.9679d0   ! +- 0.0096 using s_ij and delta_13 from PDG 2004
      rtckm2 = 0.0076d0
      lamnp2 = 0.128d0
         del = 0.900d0

C   hep-ph/9903226:

       cbue1 = 65/6.d0 - 4*zeta2
       cbue2 = (3763495.d0 - 48*(35453 + 4104*dlog(2.d0))*zeta2 
     .       - 938592*zeta3 + 936360*zeta4)/23328.d0

      call alfahat(mb,dgamma,alpem)
      alpem = alpem/pi1
      asb = alphas(mb)/pi1

C  Appendix C of hep-ph/9802391:

      a7(1) =   14/23.d0
      a7(2) =   16/23.d0
      a7(3) =    6/23.d0
      a7(4) = - 12/23.d0
      a7(5) =   0.4086d0
      a7(6) = - 0.4230d0
      a7(7) = - 0.8994d0
      a7(8) =   0.1456d0

C     h3(1) =           0.0000d0
C     h3(2) =           0.0000d0
C     h3(3) =        2/    63.d0
C     h3(4) = -      1/    27.d0
C     h3(5) = -         0.0659d0
C     h3(6) =           0.0595d0
C     h3(7) = -         0.0218d0
C     h3(8) =           0.0335d0

C     h4(1) =           0.0000d0
C     h4(2) =           0.0000d0
C     h4(3) =        1/    21.d0
C     h4(4) =        1/     9.d0
C     h4(5) =           0.0237d0
C     h4(6) = -         0.0173d0
C     h4(7) = -         0.1336d0
C     h4(8) = -         0.0316d0

C     h5(1) =           0.0000d0
C     h5(2) =           0.0000d0
C     h5(3) = -      1/   126.d0
C     h5(4) =        1/   108.d0
C     h5(5) =           0.0094d0
C     h5(6) = -         0.0100d0
C     h5(7) =           0.0010d0
C     h5(8) = -         0.0017d0

C     h6(1) =           0.0000d0
C     h6(2) =           0.0000d0
C     h6(3) = -      1/    84.d0
C     h6(4) = -      1/    36.d0
C     h6(5) =           0.0108d0
C     h6(6) =           0.0163d0
C     h6(7) =           0.0103d0
C     h6(8) =           0.0023d0

      h7(1) =   626126/272277.d0
      h7(2) = -  56281/ 51730.d0
      h7(3) = -      3/     7.d0
      h7(4) = -      1/    14.d0
      h7(5) = -         0.6494d0
      h7(6) = -         0.0380d0
      h7(7) = -         0.0186d0
      h7(8) = -         0.0057d0

      h8(1) =   313063/363036.d0
      h8(2) =           0.0000d0
      h8(3) =           0.0000d0
      h8(4) =           0.0000d0
      h8(5) = -         0.9135d0
      h8(6) =           0.0873d0
      h8(7) = -         0.0571d0
      h8(8) =           0.0209d0

      e7(1) =   4661194/816831.d0
      e7(2) = -    8516/  2217.d0
      e7(3) =                0.d0
      e7(4) =                0.d0
      e7(5) = -          1.9043d0
      e7(6) = -          0.1008d0
      e7(7) =            0.1216d0
      e7(8) =            0.0183d0

      f7(1) = - 17.3023d0
      f7(2) =    8.5027d0
      f7(3) =    4.5508d0
      f7(4) =    0.7519d0
      f7(5) =    2.0040d0
      f7(6) =    0.7476d0
      f7(7) = -  0.5385d0
      f7(8) =    0.0914d0

      kl7(1) =   9.93719d0 + 15*0.5783914d0
      kl7(2) = - 7.48777d0 - 15*0.3921387d0
      kl7(3) =   1.26884d0 - 15*0.1428571d0
      kl7(4) = - 0.29251d0 + 15*0.0476190d0
      kl7(5) = - 2.29231d0 - 15*0.1274610d0
      kl7(6) = - 0.14613d0 + 15*0.0317363d0
      kl7(7) =   0.12387d0 + 15*0.0078114d0
      kl7(8) =   0.08123d0 - 15*0.0031013d0

C   Semileptonic phase space factor according to Y. Nir, PLB 221 (1989) 184.
C   hphsp takes into account the common scale in z, and assumes the MS-bar
C   mass for the mb**5 factor multiplying the semileptonic decay. 

          z  = (mcrun(mb)/mb)**2
          z2 = z*z
          z3 = z*z2
          z4 = z*z3
          z5 = z*z4
          z6 = z*z5
       logz  = dlog(z)
       logz2 = logz**2
       logz3 = logz*logz2
         l3p = logz3 + pi2*logz
      rootz  = dsqrt(z)

      hphsp = (z2 - 1)*(65*(z2 + 1)/4 - 479*z/3)
     .      +  z*logz*( - 4 + 174*z + 212*z2/3 - 19*z3/3)
     .      + z2*logz2*(z2 - 36) - 32*rootz**3*(1 + z)
     .      * (pi2 - 2*logz*dlog((1 - rootz)/(1 + rootz))
     .      - 4*dreal(spence(rootz*(1,0)) - spence(rootz*( - 1,0))))
     .      + (1 - z2)*(17*(z2 + 1)/3 - 64*z/3)*dlog(1 - z)
     .      - 4*(1 + 30*z2 + z4)*logz*dlog(1 - z)
     .      - (1 + 16*z2 + z4)*(6*dreal(spence(z*(1,0))) - pi2)

      fphsp = 1 - 8*z + 8*z3 - z4 - 12*z2*logz - 2*asb/3*hphsp

C  r7 assumes MS-bar mass for the mb**5 factor multiplying the partial width.

        del2 = del*del
        del3 = del*del2
      logdel = dlog(del)
      logdl2 = dlog(1 - del)
      li2del = dreal(spence((1 - del)*(1,0)))
      sudfac = dexp( - asb*(2*logdel**2 + 7*logdel)/3)

      r2 = - 2*(833 - 144*pi2*rootz**3 - 36*(48 - 5*pi2 - 36*zeta3 
     .   + 9*(4 - pi2)*logz + 3*logz2 + logz3)*z
     .   - 36*(18 + 2*pi2 + 6*(2 - pi2)*logz + logz3)*z2
     .   + 6*(9 + 14*pi2 - 182*logz + 126*logz2)*z3)/243
     .   + (0,16)*((28 - 12*logz)*z3 + 3*(3*logz2 - pi2)*z2
     .   + 3*(15 - pi2 + 3*logz*(logz + 1))*z - 5)*pi1/81
      r7 = (42 - 8*pi2)/9
      r8 = - 4*(2*pi2 - 33 - (0,6)*pi1)/27

C  f221 and f271 are the (exact) integrals in Appendix E of hep-ph/9802391
C  for t = 0..4.  The remainder, t = 4..1/z, can also be calculated exactly,
C  but the result contains polylogarithms up to 5th order and arccosh
C  and is very complicated; here we use an expansion in z. 

      f221 = 8*((2 - 14*zeta3 + 4*dlog(2.d0)*pi2 - 45*zeta4/4)*z 
     .     - (24*(1 - zeta2) - 18*pi2*zeta3 + 180*dlog(2.d0)*zeta4 
     .     + 93*zeta5)*z2 + (224/3.d0 - 17*pi2 + 90*zeta4)*z3)/27
      f271 = 8*(pi2 - 8)/9*z2 + 4*(112 - 15*pi2)/27*z3
      f22  = f221 + (4 + (36*pi2 + 96*pi2*dlog(2.d0) + 270*zeta4 + 
     .       144*zeta3 + 24*l3p - 36*logz2 - 32*logz3 - 84*logz - 138)*z
     .     + 8*z2*(12*(6 + 9*zeta2*zeta3 + 6*logz + l3p + logz2 - 2*pi2) 
     .     + 135*(2*zeta4*dlog(4*z) + zeta5) + 2*logz3*(3*logz2/10+pi2))
     .     + 4*(300*logz - 565 - 18*pi2 - 1080*zeta4 - 144*logz2 -
     .       6*logz2*(logz2 + 2*pi2) + 36*l3p)*z3
     .     + (3716/9.d0 + 56*pi2 - 1304*logz/3 - 48*l3p + 168*logz2)*z4
     .     + (115/3.d0 + 44*pi2 - 664*logz - 240*l3p + 132*logz2)/9*z5
     .     - (26011/20.d0 + 330*(pi2 + 3*logz2) + 3681*logz + 
     .       1400*l3p)/50*z6)/81
      f27  = f271 + ( - 2 + 3*(2*pi2 - 7 - 6*logz - 2*logz2)*z 
     .     + 24*(8 - 2*pi2 + 2*logz + logz2)*z2 
     .     + (96*pi2 - 382 + 24*logz - 36*logz2)*z3 + 16*(3 - 5*logz)*z4 
     .     - (13 + 70*logz)*z5 - 8*(32 + 63*logz)/5*z6)/27
      f28  = - f27/3
      f77  = (10*del + del2 - 2*del3/3 + del*(del - 4)*logdel)/3
      f78  = 8*(li2del - zeta2 - del*logdel + del/4*(9-del + del2/3))/9
      f88  = (4*li2del - 2*pi2/3 + 8*logdl2 - del*(2 + del)*logdel 
     .     + 7*del + 3*del2 - 2*del3/3 - 2*(2*del + del2 
     .     + 4*logdl2)*dlog(mb/msrun(mb)))/27

       x  = mtrun(mw)**2/mw2
       x2 = x*x
       x3 = x*x2
       x4 = x*x3
       x5 = x*x4
      x12 = (x - 1)**2
      x13 = (x - 1)*x12
      x14 = (x - 1)*x13
      x15 = (x - 1)*x14
      lnx = dlog(x)
      li2 = dreal(spence((1 - 1/x)*(1,0)))

C  Two loop EW corrections from hep-ph/9804252:

      ftw = pi2*(2 - 3*x)*(2 - x)*x/48
     .    - x*(176 - 373*x + 491*x2 - 468*x3 + 72*x4)/192/x13
     .    + (x - 2)*x2*(5 - 14*x + 11*x2 - 3*x3)/8/x13*li2
     .    - x*(2 + x)*(7 + 16*x - 47*x2)/96/x12*dlog(x - 1)
     .    + x*(80 - 115*x + 200*x2 - 425*x3 + 220*x4 - 83*x5)/96/x14
     .      *lnx - x2*(4 - 5*x - 3*x2)/16/x13*lnx*
     .      ((5 - 3*x + x3)/x12*lnx - (2 + x)*dlog(x - 1))

C  Appendix A of hep-ph/9802391:

      C70 = x*( - 8*x3 + 3*x2 + 12*x - 7 + 6*x*(3*x - 2)*lnx)/24/x14
      C80 = x*( -   x3 + 6*x2 -  3*x - 2 - 6*x          *lnx)/ 8/x14
      C41 = x*(x2 + 11*x - 18)/12/x13 - 2*(lnx + 1)/3
     .    + x2*(4*x2 - 16*x + 15)/6/x14*lnx
      C71 =  x*( - 16*x3 - 122*x2 + 80*x -  8)/9/x14*li2
     .    + x2*(             6*x2 + 46*x - 28)/3/x15*lnx**2 - lnx*
     .   (102*x5 +  588*x4 +  2262*x3 -  3244*x2 + 1364*x - 208)/ 81/x15
     .    +       (1646*x4 + 12205*x3 - 10740*x2 + 2509*x - 436)/486/x14
      C81 =  x*( - 4*x3 + 40*x2 + 41*x +  1)/6/x14*li2
     .    - x2*(                  17*x + 31)/2/x15*lnx**2   - lnx*
     .  (210*x5 - 1086*x4 -  4893*x3 -  2857*x2 + 1994*x - 280)/ 216/x15
     .    +       (737*x4 - 14102*x3 - 28209*x2 +  610*x - 508)/1296/x14

C  Section 3 of hep-ph/9802391 and the author's program:

      eta  = alphas(mw)/pi1/asb

C     C10eff =    eta7(3) - eta7(4)
C     C20eff = (2*eta7(3) + eta7(4))/3
C     C30eff = 0.d0
C     C40eff = 0.d0
C     C50eff = 0.d0
C     C60eff = 0.d0
       sum70 = 0.d0
       sum80 = 0.d0
       sum71 = 0.d0
      do 10 i = 1, 8
         eta7(i) = eta**a7(i)
C         C30eff = C30eff + eta7(i)*h3(i)
C         C40eff = C40eff + eta7(i)*h4(i)
C         C50eff = C50eff + eta7(i)*h5(i)
C         C60eff = C60eff + eta7(i)*h6(i)
           sum70 =  sum70 + eta7(i)*h7(i)
           sum80 =  sum80 + eta7(i)*h8(i)
           sum71 =  sum71 + eta7(i)*(f7(i) + eta*(e7(i)*C41 + kl7(i)))
 10   continue

      do 20 i = 9, 12
         eta7(i) = eta7(i-8)/eta
 20   continue

C  We will use a different C20eff = C20eff - C10eff/6, as defined by
C  A.J. Buras, M. Misiak, M. Münz, S. Pokorski, NPB 424, 374 (1994).

      C20eff = (eta7(3) + eta7(4))/2 
      C70eff =  eta7(2)*C70 + 8*(eta7(1) - eta7(2))/3*C80 + sum70
      C80eff =  eta7(1)*C80                               + sum80
      C71eff = (eta7(2)*C71 + 8*(eta7(1) - eta7(2))/3*C81)*eta
     .       + 4*((  74416       - 1674721*eta/25)*eta7(2)
     .       -    (1791104/25.d0 -   64217*eta   )*eta7(1))*C80/14283
     .       + 37208*(eta - 1)*eta7(2)*C70/4761 + sum71

C  resummed 2-loop effects from hep-ph/9805303:

      C7em = (32*eta7(9)/75 - 40*eta7(10)/69 + 88*eta7(2)/575)*C70
     .     + ( - 32*eta7(9)/575 + 32*eta7(10)/1449 + 640*eta7(1)/1449 
     .         - 704*eta7(2)/1725)*C80 - 748*eta7(2)/8625 
     .     - 190*eta7(12)/8073 - 359*eta7(11)/3105 + 4276*eta7(4)/121095 
     .     + 350531*eta7(9)/1009125 + 2*eta7(10)/4347 
     .     - 5956*eta7(3)/15525 + 38380*eta7(1)/169533 

      Vrtuel = r2*C20eff + r7*C70eff + r8*C80eff
      Abrems = asb*(C20eff*(C80eff*f28 + C70eff*f27  + C20eff*f22)
     .       +      C80eff*(C80eff*f88 + C70eff*f78) + C70eff**2*f77)

      dabs2 = cdabs(C70eff + asb/4*(C71eff + Vrtuel))**2

C  Correction for phase space difference from small b --> u e nu component.

      fphsp  = fphsp + rtckm2*(1 + cbue1*asb + cbue2*asb**2 - fphsp)

C  Non-perturbative contributions:

      delnpg = - 9*lamnp2/2*(C70eff/mb)**2
      delnpc = -   lamnp2/9/mc**2*C70eff*C20eff
      dlnpsl =   3*lamnp2*(fphsp - 4*(1 - z)**4)/2/mb**2

C  Correction from non-peturbative b --> u e nu component:

      dlnpsl = dlnpsl - 6*rtckm2*lamnp2*(1 - (1 - z)**4)/mb**2

      sgcenu = 6*alpha/pi1*ratckm*sudfac*(dabs2 + Abrems/sudfac + delnpg 
     .       + delnpc - alpem*C70eff*(ftw/sinhat*eta7(2) + 2/asb*C7em))
     .       / (fphsp + dlnpsl)/(1 + 6*(1/eta - 1)/23*alpem/asb)

      return
      end
