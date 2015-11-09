      subroutine kappaf(kappa,fmin,fmax)

C   following G. Degrassi and A. Sirlin, NPB 352 (1991)  342,
C         and P. Gambino  and A. Sirlin, PRD  49 (1994) 1160.
C  
C   f is the fermion index

      implicit none
      integer fmin,fmax,f
      double precision gds,gy,taub,alphas,az
      complex*16 fds,fy,f1,dkapse,dkapfv,dkbtop,kappa(fmin:fmax)
      include 'common.f'

      call deltakappase(dkapse,.false.)

      az = alphas(mz)/pi1
      if (fasmt0.eqv..false.) az = 0.d0

      f1 = fds(1.d0)  
      fy = fds(ratzw2)
      gy = gds(ratzw2)

      do 100 f = 0, 3

         dkapfv = - alphat/4/pi1/sinhat*(4*coshat*gy +
     .              (1/2.d0 - sinhat*(1 - dabs(q(f))))*fy -
     .              (1 - 6*sinhat*dabs(q(f)) + 8*sinhat**2*q(f)**2)
     .              /4/coshat*f1)

         kappa(f) = 1 + dkapse + dkapfv

 100  continue

C   O(alpha alpha_s) vertex corrections for massless u, d, s, and c quarks 
C   from A. Czarnecki and J. H. K"uhn, PRL 77 (1997) 3955. The coefficient 
C   -0.87 is from the revised version of hep-ph/9608366 dated 25.11.97. 
C   It differs from the original preprint as well as from the PRL.
C   Special m_t dependent corrections of the same order to the Zbb vertex
C   were computed by R. Harlander, T. Seidensticker, and M. Steinhauser
C   in hep-ph/9712228 using the heavy top quark expansion (see bvertex.f)

      f1 = f1*(1 - az) + 4/3.d0*0.37d0*az
      fy = fy*(1 - az) + 4/3.d0*0.37d0*az
      gy = gy*(1 - az) -        0.87d0*az

      do 200 f = 4, 9

         dkapfv = - alphat/4/pi1/sinhat*(4*coshat*gy +
     .              (1/2.d0 - sinhat*(1 - dabs(q(f))))*fy -
     .              (1 - 6*sinhat*dabs(q(f)) + 8*sinhat**2*q(f)**2)
     .              /4/coshat*f1)

         kappa(f) = 1 + dkapse + dkapfv

 200  continue

      call bvertex(dkbtop,taub)
      kappa(9) = kappa(9) + dkbtop
      if (fa2mt4.eqv..true.) kappa(9) = kappa(9)/(1 + taub)
      if (fobliq.eqv..true.) kappa(9) = kappa(9) + Bkappa 
      
      return
      end
      
      
      subroutine deltakappase(dkapse,feonly)

C   following G. Degrassi and A. Sirlin, NPB 352 (1991)  342,
C         and P. Gambino  and A. Sirlin, PRD  49 (1994) 1160.

      implicit none
      integer f
      logical feonly
      double precision lambda,alphas,al,az,at,c(2,3),t1,t2,t3
      double precision mbrun,mcrun,rtauz2,ratcz2,ratbz2,rcs
      complex*16 aboson,imaqed,logc2,dkapse
      complex*16 const1,const2,const3,xa2mt0
      include 'common.f'

      c(1,1) =       1/      5.d0
      c(1,2) =       3/    140.d0
      c(1,3) =       1/    315.d0
      c(2,1) =     194/    405.d0
      c(2,2) =     983/  18900.d0
      c(2,3) =   12079/1984500.d0

          al = alphat/pi1
          az = alphas(mz)/pi1
          at = alphas(mt)/pi1
      const1 = (5/3.d0 + (0,1)*pi1)/3
      const2 = (55/12.d0 - 4*zeta3 + (0,1)*pi1)/3
      const3 = (16.d0 + (0,12)*pi1)/3
      rtauz2 = (mtau/mz)**2
      ratcz2 = (mcrun(mz)/mz)**2
      ratbz2 = (mbrun(mz)/mz)**2
         rcs = dsqrt(coshat/sinhat)

      xa2mt0 = 1.d0
      if (fa2mt0.eqv..false.) xa2mt0 = 0.d0
      if (fasmt0.eqv..false.) then
         az = 0.d0
         at = 0.d0
      endif

      dkapse = 0.d0
      do 100 f = 1, 9
         if (f.ne.6) then
            dkapse = dkapse + rcs*al*nc(f)*q(f)*v(f)*
     .           (const1 + 3/4.d0*al*q(f)**2*const2*xa2mt0)
            if (f.ge.4) dkapse = dkapse + 3*rcs*al*q(f)*v(f)*az*const2
            if (flagmf.eqv..true.) then
               if (f.eq.3) dkapse = dkapse + rcs*al*q(f)*v(f)*2*rtauz2
               if (f.eq.5) dkapse = dkapse + rcs*al*3*q(f)*v(f)*ratcz2*
     .              (2 + az*const3)
               if (f.eq.9) dkapse = dkapse + rcs*al*3*q(f)*v(f)*ratbz2*
     .              (2 + az*const3)
            endif
         else
            t1 = 1.d0/rattz2
            t2 = 1.d0/rattz2**2
            t3 = 1.d0/rattz2**3
            dkapse = dkapse + rcs*al*q(f)*v(f)*
     .           ((c(1,1)*t1 + c(1,2)*t2 + c(1,3)*t3) + at*
     .            (c(2,1)*t1 + c(2,2)*t2 + c(2,3)*t3))
          if ((flagmr.eqv..false.).or.(mt.le.mz)) 
     .           dkapse = dkapse - rcs*al*q(f)*v(f)*(dlog(rattz2)*
     .           (1 + at + al*xa2mt0/3) - 13/12.d0*at - 5*al*xa2mt0/4)
         endif
 100  continue
      
      imaqed = 0.d0
      aboson = 0.d0
       logc2 = 0.d0
      if (feonly.eqv..false.) then
         if (fla2im.eqv..true.) then
            imaqed  = alphat*20/9*dimag(dkapse)
         endif
         
         aboson = al/4/sinhat*( - (3*coshat + 1/6.d0)*dlog(ratzw2)
     .          - 1/9.d0 + (lambda(ratzw2) - 1)*
     .            (8/ratzw2*(coshat + 1/3.d0) + 6*coshat + 1/3.d0))

         logc2   = - al/sinhat*coshat*dlog(ratzw2)
      endif

      dkapse = dkapse + aboson + imaqed + logc2

      return
      end
