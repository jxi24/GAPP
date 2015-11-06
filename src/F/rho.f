      subroutine rhof(rho,fmin,fmax)
C
C   following G. Degrassi and A. Sirlin, NPB 352 (1991)  342,
C  
C   f is the fermion index
C
      implicit none
      integer fmin,fmax,f
      double precision gds,gy,hds,taub,lambda,droseb,alphas,az
      complex*16 fds,fy,f1,drhofv,rho(fmin:fmax),dkbtop,drosef
      include 'common.f'

      droseb = alphat/4/pi1/sinhat*(4*coshat*dlog(ratzw2) +
     .         5/36.d0*(4 - 1/coshat + 12/coshat/ratzw2) +
     .         3*coshat - 8/3.d0/ratzw2 - 4*coshat/ratzw2 -
     .         (1/12.d0/coshat - 1/3.d0 - 3*coshat)*dlog(ratzw2) +
     .         lambda(ratzw2)*(1/6.d0/coshat + 1/3.d0/coshat/ratzw2 -
     .         2/3.d0 - 4/3.d0/ratzw2 - 6*coshat*(1 - 1/ratzw2) +
     .         (18*coshat - 8*sinhat**2/coshat/ratzw2)/(4 - ratzw2)) +
     .         hds(rathz2)/2/coshat)

      call deltarhosef(drosef)

      az = alphas(mz)/pi1
      if (fasmt0.eqv..false.) az = 0.d0

      f1 = fds(1.d0)  
      fy = fds(ratzw2)
      gy = gds(ratzw2)

      do 100 f = 0, 3
         
         drhofv = alphat/4/pi1/sinhat*(8*coshat*gy +
     .            (1 - 2*sinhat*(1 - dabs(q(f))))*fy -
     .            (1 - 6*sinhat*dabs(q(f)) + 12*(sinhat*q(f))**2)
     .            /2/coshat*f1)
         
         rho(f) = (1 + droseb + drhofv)/(1 - drosef)

 100  continue
      
C   O(alpha alpha_s) vertex corrections for massless u, d, s, and c quarks 
C   from A. Czarnecki and J. H. K"uhn, PRL 77 (1997) 3955. The coefficient 
C   -0.87 is from the revised version of hep-ph/9608366 dated 25.11.97. 
C   It differs from the original preprint as well as from the PRL.
C   Special m_t dependent corrections of the same order to the Zbb vertex
C   were computed by R. Harlander, T. Seidensticker, and M. Steinhauser
C   in hep-ph/9712228 using the heavy top quark expansion.

      f1 = f1*(1 - az) + 4/3.d0*0.37d0*az
      fy = fy*(1 - az) + 4/3.d0*0.37d0*az
      gy = gy*(1 - az) -        0.87d0*az

      call bvertex(dkbtop,taub)

      do 200 f = 4, 9
         
         drhofv = alphat/4/pi1/sinhat*(8*coshat*gy +
     .            (1 - 2*sinhat*(1 - dabs(q(f))))*fy -
     .            (1 - 6*sinhat*dabs(q(f)) + 12*(sinhat*q(f))**2)
     .            /2/coshat*f1)
         
         if (f.ne.9) then
            rho(f) = (1 + droseb + drhofv)/(1 - drosef)
         else
            rho(9) = (1 + droseb + drhofv - 2*dkbtop)/(1 - drosef)
            if (fa2mt4.eqv..true.) rho(9) = rho(9)*(1 + taub)**2
            if (fobliq.eqv..true.) rho(9) = rho(9) + Brho
         endif
         
 200  continue
      
      return
      end


      subroutine deltarhosef(drosef)

C   following G. Degrassi and A. Sirlin, NPB 352 (1991)  342,

      implicit none
      integer f
      double precision alphas,al,az,at,c(4,3),t1,t2,t3,dfks,dfksh
      double precision mbrun,mcrun,mtauz2,mcz2,mbz2,xa2mt0
      complex*16 const1,const2,const3,const4,const5,drosef
      include 'common.f'

      c(1,1) =       2/      5.d0
      c(1,2) =       9/    140.d0
      c(1,3) =       4/    315.d0
      c(2,1) =       1/      5.d0
      c(2,2) =       3/    140.d0
      c(2,3) =       1/    315.d0
      c(3,1) =     388/    405.d0
      c(3,2) =     983/   6300.d0
      c(3,3) =   12079/ 496125.d0
      c(4,1) =     257/    810.d0
      c(4,2) =     251/  12600.d0
      c(4,3) = -  1187/1984500.d0

          al = alphat/pi1
          az = alphas(mz)/pi1
          at = alphas(mt)/pi1
      const1 = - 2/9.d0
      const2 = - (43/12.d0 - 4*zeta3)/3.d0 
      const3 = -       2.d0*(1.d0 + (0,1)*pi1)
      const4 =         4.d0*(1.d0 + (0,1)*pi1)
      const5 = - 22.d0/3.d0*(1.d0 + (0,1)*pi1) - 4.d0*(0,1)*pi1
      mtauz2 = (mtau/mz)**2
        mcz2 = (mcrun(mz)/mz)**2
        mbz2 = (mbrun(mz)/mz)**2

      xa2mt0 = 1.d0
      if (fa2mt0.eqv..false.) xa2mt0 = 0.d0
      if (fasmt0.eqv..false.) then
         az = 0.d0
         at = 0.d0
      endif

      drosef = 0.d0
      do 100 f = 0, 9
         if (f.ne.6.) then
            drosef = drosef + al*nc(f)*(v(f)**2 + a(f)**2)*
     .           (const1 + 3/4.d0*al*q(f)**2*const2*xa2mt0)
            if (f.ge.4) drosef = drosef + 
     .           al*3.d0*(v(f)**2 + a(f)**2)*az*const2
            if (flagmf.eqv..true.) then
               if (f.eq.3) drosef = drosef + al*a(f)**2*mtauz2*const3
               if (f.eq.5) drosef = drosef + al*mcz2*3.d0*
     .              (v(f)**2*az*const4 + a(f)**2*(const3 + az*const5))
               if (f.eq.9) drosef = drosef + al*mbz2*3.d0*
     .              (v(f)**2*az*const4 + a(f)**2*(const3 + az*const5))
            endif
         else
            t1 = 1.d0/rattz2
            t2 = 1.d0/rattz2**2
            t3 = 1.d0/rattz2**3
            drosef = drosef - al*v(f)**2*
     .           ((c(1,1)*t1 + c(1,2)*t2 + c(1,3)*t3) + at*
     .            (c(3,1)*t1 + c(3,2)*t2 + c(3,3)*t3)) - al*a(f)**2*
     .           ((c(2,1)*t1 + c(2,2)*t2 + c(2,3)*t3 - 1) + at*
     .            (c(4,1)*t1 + c(4,2)*t2 + c(4,3)*t3) -
     .            (at + al*xa2mt0/3)*17/9.d0)
              dfks = (1/sinhat/6 - 4/9.d0)*(dlog(rattz2)*(1 + at
     .             + xa2mt0*al/3) - 13*at/12 - xa2mt0*5/4*al)
             dfksh = dfks/(1 - 8*sinhat/3.d0)
            drosef = drosef + al*(v(f)**2 + a(f)**2)*6*sinhat*dfksh
            if (flagmr.eqv..true.) then
               drosef = drosef + al*(sinhat/coshat*dfks - dfksh)
            endif
         endif
 100  continue

      if (fobliq.eqv..true.) drosef = drosef 
     .                              + alphat/4/sinhat/coshat*Spar

      return
      end


      double precision function hds(x)
      
      implicit none
      double precision x,xx,omegds

      xx = x**2

      hds = 31/18.d0 - x + xx/3.d0 +
     .      (1 - 3/2.d0*x + 3/4.d0*xx - x**3/6.d0)*dlog(x) + 
     .      (2/(4 - x) - 2 + 5/6.d0*x -   xx/3.d0)*omegds(x)

      return
      end


      double precision function omegds(x)
      
      implicit none
      double precision x,root

      if (x.ge.4.d0) then
         root  = dsqrt(1 - 4/x)
         omegds = x/2.d0*root*dlog((1 - root)/(1 + root))
      else
         root  = dsqrt(4/x - 1)
         omegds = x*root*datan(root)
      endif

      return
      end
