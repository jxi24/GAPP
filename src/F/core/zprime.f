      subroutine zprimecoup(xval,npar)

      implicit none
      integer f,npar
      double precision norm,alpzLR,xval(npar)
      double precision aparam,bparam,cotalp,ratbsa
      include 'common.f'

C for the Z_chi defined through SO(10) -->  SU(5) X U(1)_chi:       flagzp = 1;
C for the Z_psi defined through    E_6 --> SO(10) X U(1)_psi:       flagzp = 2;
C for the Z_eta defined through sqrt(3/8) Z_chi - sqrt(5/8) Z_psi:  flagzp = 3;
C for the Z_LR defined through SU(2)_R X U(1) --> U(1)_Y X U(1)_LR: flagzp = 4;
C for the Z^prime defined to have the same couplings as the Z:      flagzp = 5;
C for the Z^string:                                                 flagzp = 6;
C for the general E_6 boson:                                        flagzp = 7;
C for the general family-universal Z^prime boson:                   flagzp = 8;
C for the model independent Z^prime boson:                          flagzp = 9;
C for the model independent Z^prime boson in v,a parametrization:   flagzp =10;
C for the general E_6 boson in alpha, beta parametrization:         flagzp =11;
C for the twisted E_6 boson:                                        flagzp =12.
C for the seesaw neutrino mass related one (Adhikari, Erler, Ma):   flagzp =13.

      if (flagzp.eq.1) then

         norm = 2*dsqrt(10.d0)
         
         eps2_L(0) =  3/norm
         eps2_L(1) =  3/norm
         eps2_L(2) =  3/norm
         eps2_L(3) =  3/norm
         eps2_L(4) = -1/norm
         eps2_L(5) = -1/norm
         eps2_L(6) = -1/norm
         eps2_L(7) = -1/norm
         eps2_L(8) = -1/norm
         eps2_L(9) = -1/norm
         
         eps2_R(0) =  5/norm
         eps2_R(1) =  1/norm
         eps2_R(2) =  1/norm
         eps2_R(3) =  1/norm
         eps2_R(4) =  1/norm
         eps2_R(5) =  1/norm
         eps2_R(6) =  1/norm
         eps2_R(7) = -3/norm
         eps2_R(8) = -3/norm
         eps2_R(9) = -3/norm

      else if (flagzp.eq.2) then
         
         eps2_L(0) =  1/dsqrt(24.d0)
         eps2_R(0) = -1/dsqrt(24.d0)
         do 10 f = 1, 9
            eps2_L(f) = eps2_L(0)
            eps2_R(f) = eps2_R(0)
 10      continue

      else if (flagzp.eq.3) then

         norm = 2*dsqrt(15.d0)
         
         eps2_L(0) =  1/norm
         eps2_L(1) =  1/norm
         eps2_L(2) =  1/norm
         eps2_L(3) =  1/norm
         eps2_L(4) = -2/norm
         eps2_L(5) = -2/norm
         eps2_L(6) = -2/norm
         eps2_L(7) = -2/norm
         eps2_L(8) = -2/norm
         eps2_L(9) = -2/norm
         
         eps2_R(0) =  5/norm
         eps2_R(1) =  2/norm
         eps2_R(2) =  2/norm
         eps2_R(3) =  2/norm
         eps2_R(4) =  2/norm
         eps2_R(5) =  2/norm
         eps2_R(6) =  2/norm
         eps2_R(7) = -1/norm
         eps2_R(8) = -1/norm
         eps2_R(9) = -1/norm

      else if (flagzp.eq.4) then

         norm = 2*dsqrt(5/3.d0)
         alpzLR = dsqrt(ratgRL**2*coshat/sinhat - 1)

         eps2_L(0) =  1/norm/alpzLR
         eps2_L(1) =  1/norm/alpzLR
         eps2_L(2) =  1/norm/alpzLR
         eps2_L(3) =  1/norm/alpzLR
         eps2_L(4) = -1/norm/alpzLR/3
         eps2_L(5) = -1/norm/alpzLR/3
         eps2_L(6) = -1/norm/alpzLR/3
         eps2_L(7) = -1/norm/alpzLR/3
         eps2_L(8) = -1/norm/alpzLR/3
         eps2_L(9) = -1/norm/alpzLR/3
         
         do 20 f = 0, 9
            eps2_R(f) = eps2_L(f) + 2*i3(f)*alpzLR/norm
 20      continue

      else if (flagzp.eq.5) then

         do 30 f = 0, 9
            eps2_L(f) =  i3(f) - q(f)*sinhat
            eps2_R(f) =        - q(f)*sinhat
 30      continue

      else if (flagzp.eq.6) then

         norm = 100.d0
         
         eps2_L(0) = - 65/norm
         eps2_L(1) = -204/norm
         eps2_L(2) = - 65/norm
         eps2_L(3) =   74/norm
         eps2_L(4) =   68/norm
         eps2_L(5) =   68/norm
         eps2_L(6) =  -71/norm
         eps2_L(7) =   68/norm
         eps2_L(8) =   68/norm
         eps2_L(9) =  -71/norm
         
         eps2_R(0) =    0/norm
         eps2_R(1) =    9/norm
         eps2_R(2) =    9/norm
         eps2_R(3) = -130/norm
         eps2_R(4) = -  6/norm
         eps2_R(5) = -  6/norm
         eps2_R(6) =  133/norm
         eps2_R(7) =    3/norm
         eps2_R(8) =    3/norm
         eps2_R(9) = -136/norm

      else if (flagzp.eq.7) then

         norm = dsqrt(16*xval(19)*xval(20)/3 + 20*xval(19)*xval(21)/3
     .        +                                 4*xval(20)*xval(21)
     .        + 40*xval(19)**2/9 + 4*xval(20)**2 + 16*xval(21)**2)

         eps2_L(1) = (-   xval(19)                          )/norm
         eps2_L(2) = (-   xval(19)                          )/norm
         eps2_L(3) = (-   xval(19)                          )/norm
         eps2_L(4) = (    xval(19)/3            +   xval(21))/norm
         eps2_L(5) = (    xval(19)/3            +   xval(21))/norm
         eps2_L(6) = (    xval(19)/3            +   xval(21))/norm
         eps2_L(7) = (    xval(19)/3            +   xval(21))/norm
         eps2_L(8) = (    xval(19)/3            +   xval(21))/norm
         eps2_L(9) = (    xval(19)/3            +   xval(21))/norm
         eps2_L(0) = (-   xval(19)                          )/norm
         
         eps2_R(1) = (-   xval(19)/3 + xval(20) -   xval(21))/norm
         eps2_R(2) = (-   xval(19)/3 + xval(20) -   xval(21))/norm
         eps2_R(3) = (-   xval(19)/3 + xval(20) -   xval(21))/norm
         eps2_R(4) = (-   xval(19)/3 - xval(20) -   xval(21))/norm
         eps2_R(5) = (-   xval(19)/3 - xval(20) -   xval(21))/norm
         eps2_R(6) = (-   xval(19)/3 - xval(20) -   xval(21))/norm
         eps2_R(7) = (    xval(19)   + xval(20)             )/norm
         eps2_R(8) = (    xval(19)   + xval(20)             )/norm
         eps2_R(9) = (    xval(19)   + xval(20)             )/norm
         eps2_R(0) = (- 5*xval(19)/3 - xval(20) - 2*xval(21))/norm

C         eps2_L(3) = (+ 2*xval(19)/3 + xval(20) -   xval(21))/norm
C         eps2_R(9) = (- 2*xval(19)/3            +   xval(21))/norm

      else if ((flagzp.eq.11).or.(flagzp.eq.12)) then

         cotalp = 1/tan(xval(19))
         ratbsa = tan(xval(20))/sin(xval(19))

         if (xval(19).eq.0.d0) then
            aparam = 0.d0
            bparam = -4*tan(xval(20))/(tan(xval(20)) + dsqrt(27/5.d0))/3
         else
            aparam = 10/(dsqrt(10.d0)*ratbsa + dsqrt(54.d0)*cotalp - 6)
            bparam = - dsqrt(8/45.d0)*aparam*ratbsa
         endif

         norm = dsqrt(16*aparam/3 + 20*bparam/3 + 4*aparam*bparam
     .        + 40/9.d0 + 4*aparam**2 + 16*bparam**2)

         eps2_L(1) = (- 1                         )/norm
         eps2_L(2) = (- 1                         )/norm
         eps2_L(3) = (- 1                         )/norm
         eps2_L(4) = (  1/3.d0          +   bparam)/norm
         eps2_L(5) = (  1/3.d0          +   bparam)/norm
         eps2_L(6) = (  1/3.d0          +   bparam)/norm
         eps2_L(7) = (  1/3.d0          +   bparam)/norm
         eps2_L(8) = (  1/3.d0          +   bparam)/norm
         eps2_L(9) = (  1/3.d0          +   bparam)/norm
         eps2_L(0) = (- 1                         )/norm
         
         eps2_R(1) = (- 1/3.d0 + aparam -   bparam)/norm
         eps2_R(2) = (- 1/3.d0 + aparam -   bparam)/norm
         eps2_R(3) = (- 1/3.d0 + aparam -   bparam)/norm
         eps2_R(4) = (- 1/3.d0 - aparam -   bparam)/norm
         eps2_R(5) = (- 1/3.d0 - aparam -   bparam)/norm
         eps2_R(6) = (- 1/3.d0 - aparam -   bparam)/norm
         eps2_R(7) = (  1      + aparam           )/norm
         eps2_R(8) = (  1      + aparam           )/norm
         eps2_R(9) = (  1      + aparam           )/norm
         eps2_R(0) = (- 5/3.d0 - aparam - 2*bparam)/norm

         if (flagzp.eq.12) then
            eps2_R(7) = (- 2/3.d0            + bparam)/norm
            eps2_R(8) = (- 2/3.d0            + bparam)/norm
            eps2_R(9) = (- 2/3.d0            + bparam)/norm
         endif

      else if (flagzp.eq.8) then

         eps2_L(1) = xval(19)
         eps2_L(2) = xval(19)
         eps2_L(3) = xval(19)
         eps2_L(4) = xval(21)
         eps2_L(5) = xval(21)
         eps2_L(6) = xval(21)
         eps2_L(7) = xval(21)
         eps2_L(8) = xval(21)
         eps2_L(9) = xval(21)
         eps2_L(0) = xval(19)
         
         eps2_R(1) = xval(20)
         eps2_R(2) = xval(20)
         eps2_R(3) = xval(20)
         eps2_R(4) = xval(22)
         eps2_R(5) = xval(22)
         eps2_R(6) = xval(22)
         eps2_R(7) = xval(23)
         eps2_R(8) = xval(23)
         eps2_R(9) = xval(23)
         eps2_R(0) = 0.d0

      else if (flagzp.eq.9) then

         eps2_L(1) = xval(19)
         eps2_L(2) = xval(19)
         eps2_L(3) = xval(24)
         eps2_L(4) = xval(21)
         eps2_L(5) = xval(21)
         eps2_L(6) = xval(26)
         eps2_L(7) = xval(21)
         eps2_L(8) = xval(21)
         eps2_L(9) = xval(26)
         eps2_L(0) = (2*xval(19) + xval(24))/3
         
         eps2_R(1) = xval(20)
         eps2_R(2) = xval(20)
         eps2_R(3) = xval(25)
         eps2_R(4) = xval(22)
         eps2_R(5) = xval(22)
         eps2_R(6) = 0.d0
         eps2_R(7) = xval(23)
         eps2_R(8) = xval(23)
         eps2_R(9) = xval(27)
         eps2_R(0) = 0.d0

      endif

      do 40 f = 0, 9 
         v2(f) = eps2_L(f) + eps2_R(f)
         a2(f) = eps2_L(f) - eps2_R(f)
 40   continue

      if (flagzp.eq.10) then

         v2(1) = xval(19)
         v2(2) = xval(19)
         v2(3) = xval(24)
         v2(4) = xval(21)
         v2(5) = xval(21)
         v2(6) = 0.d0
         v2(7) = xval(22)
         v2(8) = xval(22)
         v2(9) = xval(26)
         v2(0) = (xval(19) + xval(20) + (xval(24) + xval(25))/2)/3

         a2(1) = xval(20)
         a2(2) = xval(20)
         a2(3) = xval(25)
         a2(4) = xval(23)
         a2(5) = xval(23)
         a2(6) = 0.d0
         a2(7) = xval(21) - xval(22) + xval(23)
         a2(8) = xval(21) - xval(22) + xval(23)
         a2(9) = xval(27)
         a2(0) = (xval(19) + xval(20) + (xval(24) + xval(25))/2)/3

         do 50 f = 0, 9 
            eps2_L(f) = (v2(f) + a2(f))/2
            eps2_R(f) = (v2(f) - a2(f))/2
 50      continue

      endif

      if (flagzp.eq.13) then

         norm = 1.d0
         
         eps2_L(0) =  xval(19)/norm
         eps2_L(1) =  xval(19)/norm
         eps2_L(2) =  xval(19)/norm
         eps2_L(3) =  xval(19)/norm
         eps2_L(4) =  xval(20)/norm
         eps2_L(5) =  xval(20)/norm
         eps2_L(6) =  xval(20)/norm
         eps2_L(7) =  xval(20)/norm
         eps2_L(8) =  xval(20)/norm
         eps2_L(9) =  xval(20)/norm
         
C        eps2_R(0) =  (3*xval(19) + 9*xval(20))/4/norm  ! model (B)
         eps2_R(0) =  (3*xval(19) + 9*xval(20))/8/norm  ! model (C)
         eps2_R(1) =  (5*xval(19) - 9*xval(20))/4/norm
         eps2_R(2) =  (5*xval(19) - 9*xval(20))/4/norm
         eps2_R(3) =  (5*xval(19) - 9*xval(20))/4/norm
         eps2_R(4) = -(3*xval(19) - 7*xval(20))/4/norm
         eps2_R(5) = -(3*xval(19) - 7*xval(20))/4/norm
         eps2_R(6) = -(3*xval(19) - 7*xval(20))/4/norm
         eps2_R(7) =  (3*xval(19) +   xval(20))/4/norm
         eps2_R(8) =  (3*xval(19) +   xval(20))/4/norm
         eps2_R(9) =  (3*xval(19) +   xval(20))/4/norm

      do 60 f = 0, 9 
         v2(f) = eps2_L(f) + eps2_R(f)
         a2(f) = eps2_L(f) - eps2_R(f)
 60   continue

         endif
      return
      end
