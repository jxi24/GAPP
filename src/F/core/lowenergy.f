      subroutine rho0

C   following W. J. Marciano and A. Sirlin, PRD 22 (1980) 2695.

      implicit none
      double precision alh,at,alphas,xt,cqcd1,cqcd2
      double precision rtauz2,ratcz2,ratbz2,ratbt2,mbrun
      include 'common.f'

      alh =     alphat/pi1
      at  = alphas(mt)/pi1
      xt = gf*mt2/2/pi2/dsqrt(2.d0)

      rtauz2 =      (mtau/mz)**2
      ratcz2 =        (mc/mz)**2
      ratbz2 = (mbrun(mt)/mz)**2
      ratbt2 = ratbz2/rattz2

      cqcd1  = zeta2 - 3/2.d0
      cqcd2  = 2.9772d0

      if (mh.ne.mz) then
         rhonc = 1 + alh/16/sinhat*(3*rattz2/coshat - 7 
     .        - 3*dlog(ratzw2)/sinhat + 3*rathz2/coshat
     .        * (dlog(rathz2)/(1 - rathz2) - dlog(rathw2)/(1 - rathw2)))
      else
         rhonc = 1 + alh/16/sinhat*(3*rattz2/coshat - 7
     .        - 3*dlog(ratzw2)/sinhat - 3*rathz2/coshat
     .        * (1 + dlog(rathw2)/(1 - rathw2)))
      endif

      if ((fasmt2.eqv..true.).and.(ffermi.eqv..false.)) rhonc = rhonc
     .     - alh*at/4/sinhat/coshat*rattz2*cqcd1

      if ((fas2mt.eqv..true.).and.(ffermi.eqv..false.)) rhonc = rhonc 
     .     - alh*at**2/4/sinhat/coshat*rattz2*cqcd2

      if (ffermi.eqv..false.) rhonc = rhonc 
     .     + 3*(alh/16/sinhat/coshat*rattz2)**2*rho2

      if (ffermi.eqv..true.) then
         if (fa2mt2.eqv..false.) then
            rhonc = rhonc - 3*alh/16/sinhat/coshat*rattz2 
     .                    + 3*xt/4*rhohat*(1 + xt*rhohat/4*rho2)
         else
            rhonc = rhonc + 3*(xt*rhohat)**2/16*rho2
         endif
         if (fasmt2.eqv..true.) rhonc = rhonc -    at*xt*rhohat*cqcd1
         if (fas2mt.eqv..true.) rhonc = rhonc - at**2*xt*rhohat*cqcd2
      endif

      if (fobliq.eqv..true.) rhonc = rhonc/(1 - alphat*Tpar)

      if (flagmf.eqv..true.) rhonc = rhonc + alh/16/sinhat/coshat
     .     *(3*ratbz2*(1 + 2*dlog(ratbt2)) + 3*ratcz2 + rtauz2)

      return
      end


      subroutine sin2theta0

C   Includes process independent corrections for low energy NC observables.
C   Our definition is the sin2th counterpart to the fine structure constant
C   and is similar to the one by A. Sirlin, NPB 332 (1990) 20.

      implicit none
      integer i
      double precision al,alh,az,at,alphas,xa2mt0,dkhad5,dahad5,sin2t0
      double precision cqcd,dgama5,lamda2,x,deltau,deltas,deltac,deltab
      include 'common.f'

      al  =  alpha/pi1
      alh = alphat/pi1
      az  = alphas(mz)/pi1
      at  = alphas(mt)/pi1

      xa2mt0 = 1.d0
      if (fa2mt0.eqv..false.) xa2mt0 = 0.d0
      if (fasmt0.eqv..false.) then
         az = 0.d0
         at = 0.d0
      endif

      cqcd = 55/12.d0 - 4*zeta3

      dahad5 = 0.0278d0
      dgama5 = dahad5/al + 11/9.d0*(5/3.d0 + az*cqcd)      
           x = (dgama5 - 143*az/108)/(1 + az) 

      deltab =   dlog(mz2/mb**2)/9 
      deltac = 4*dlog(mz2/mc**2)/9
      deltas =   dlog(mz2/ms**2)/9

      do 100, i=1,2
         lamda2 = mz2*dexp( - 9/5.d0*(x - deltab - deltac - deltas))
         deltas = dlog(mz2/(ms**2 + lamda2))/9
 100  continue

      deltau = 4/5.d0*(x - deltab - deltac - deltas)

      dkhad5 = al*((3/4.d0 - sinhat)*dgama5 - 13*az/36
     .       - 3*(1 + az)*(deltau + deltac)/8)

      sin2t0 = sinhat + dkhad5 + al*((1 - 4*sinhat)/12
     .       * ((dlog(mz2/me**2) + dlog(mz2/mmu**2) + dlog(mz2/mtau**2))
     .       * (1 + 3/4.d0*alh*xa2mt0) + 135/32.d0*alh*xa2mt0)
     .       - (7*coshat/4 + 1/24.d0)*dlog(ratzw2) + sinhat/6 - 7/18.d0)

      if ((flagmr.eqv..false.).or.(mt.le.mz))
     .     sin2t0 = sin2t0 - al/6*(1 - 8*sinhat/3)*(dlog(rattz2)
     .     *(1 + at + alh*xa2mt0/3) - 13/12.d0*at - 5*alh*xa2mt0/4)

      return 
      end


      double precision function shat(x)

      implicit none
      double precision x,dgamma,mudbar,msbar,mcbar,al,alh
      double precision alfax,alfamt,alfamw,alfamb,alftau,alfamc,amcbar
      double precision amsbar,audbar,alfamu,alfame,atmins,awmins,abmins
      double precision ataumi,acmins,amumin,aemins,dah3mc,dahad2,xa2mt0
      double precision alphas,ax2,ac2,atau2,ab2,aw2,az2,at,at2,awpi
      include 'common.f'

      mudbar = 0.176d0
      msbar  = 0.305d0
      mcbar  = 1.176d0

      dah3mc = dahad3*1.175d0
      dahad2 = dahad3*0.298d0

      al   = alpha/pi1
      alh  = alphat/pi1

      ax2   = alphas(x)**2/pi2
      ac2   = alphas(mc)**2/pi2
      atau2 = alphas(mtau)**2/pi2
      ab2   = alphas(mb)**2/pi2
      aw2   = alphas(mw)**2/pi2
      az2   = alphas(mz)**2/pi2
      at    = alphas(mt)/pi1
      at2   = at*at

      xa2mt0 = 1.d0
      if (fa2mt0.eqv..false.) xa2mt0 = 0.d0
      if (fasmt0.eqv..false.) then
         at = 0.d0
      endif

      call alfahat(x,dgamma,alfax)
      if (x.ge.mw) then
         if (x.lt.mt.and.flagmr.eqv..true.) then
            shat = alfax/alphat*sinhat + 21*(1 - alfax/alphat)/44
     .           + alfax/pi1*(625*dlog(x/mz)/264 
     .           + 9*dlog(alfax/alphat)/22
     .           + 5*(11 - 24*zeta3)/6072*(az2 - ax2))
         else if (flagmr.eqv..false.) then
            shat = alfax/alphat*sinhat + 9*(1 - alfax/alphat)/20
     .           + alfax/pi1*(289*dlog(x/mz)/120
     .           + 21*dlog(alfax/alphat)/110
     .           + (11 - 24*zeta3)/336*(az2 - ax2))         
         else
            call alfahat(mt - 1.d-12,dgamma,atmins)
            shat = atmins/alphat*sinhat + 21*(1 - atmins/alphat)/44
     .           + atmins/pi1*(625*dlog(mt/mz)/264 
     .           + 9*dlog(atmins/alphat)/22
     .           + 5*(11 - 24*zeta3)/6072*(az2 - at2))
            call alfahat(mt,dgamma,alfamt)
            shat = alfamt/atmins*shat + 3*(1 - alfamt/atmins)/8
            shat = alfax/alfamt*shat + 9*(1 - alfax/alfamt)/20
     .           + alfax/pi1*(289*dlog(x/mt)/120
     .           + 21*dlog(alfax/alfamt)/110
     .           + (11 - 24*zeta3)/336*(at2 - ax2))
         endif
      else
         call alfahat(mw,dgamma,alfamw)
         call alfahat(mw - 1.d-14,dgamma,awmins)
         if (flagmr.eqv..false.) then
            shat = alfamw/alphat*sinhat + 9*(1 - alfamw/alphat)/20
     .           + alfamw/pi1*(289*dlog(mw/mz)/120 
     .           + 21*dlog(alfamw/alphat)/110
     .           + (11 - 24*zeta3)/336*(az2 - aw2))
            awpi = awmins/pi1
            shat = shat - awpi/6*(1 - 8*shat/3)*(dlog(rattw2)
     .           *(1 + at + awpi*xa2mt0/3) - 13*at/12 - 5*awpi*xa2mt0/4)
     .           - awpi/6*(1 - shat)
         else
            shat = alfamw/alphat*sinhat + 21*(1 - alfamw/alphat)/44
     .           + alfamw/pi1*(625*dlog(mw/mz)/264 
     .           + 9*dlog(alfamw/alphat)/22
     .           + 5*(11 - 24*zeta3)/6072*(az2 - aw2))
            call alfahat(mw - 1.d-14,dgamma,awmins)
            shat = 1 - awmins/alfamw*(1 - shat)
         endif
         if (x.ge.mb) then
            shat = alfax/awmins*shat + 21*(1 - alfax/awmins)/44
     .           + alfax/pi1*(15*dlog(x/mw)/33 
     .           + 153*dlog(alfax/awmins)/1760
     .           + 5*(11 - 24*zeta3)/6072*(aw2 - ax2))
         else
            call alfahat(mb,dgamma,alfamb)
            shat = alfamb/awmins*shat + 21*(1 - alfamb/awmins)/44
     .           + alfamb/pi1*(15*dlog(mb/mw)/33 
     .           + 153*dlog(alfamb/awmins)/1760
     .           + 5*(11 - 24*zeta3)/6072*(aw2 - ab2))
            call alfahat(mb - 1.d-14,dgamma,abmins)
            shat = abmins/alfamb*shat + 3*(1 - abmins/alfamb)/4
            if (x.ge.mtau) then
               shat = alfax/abmins*shat + 9*(1 - alfax/abmins)/20
     .              + alfax/pi1*(2*dlog(x/mb)/5 
     .              + 3*dlog(alfax/abmins)/38
     .              + (11 - 24*zeta3)/900*(ab2 - ax2))
            else
               call alfahat(mtau,dgamma,alftau)
               shat = alftau/abmins*shat + 9*(1 - alftau/abmins)/20
     .              + alftau/pi1*(2*dlog(mtau/mb)/5 
     .              + 3*dlog(alftau/abmins)/38
     .              + (11 - 24*zeta3)/900*(ab2 - atau2))
               call alfahat(mtau - 1.d-14,dgamma,ataumi)
               shat = ataumi/alftau*shat + (1 - ataumi/alftau)/4
               if (x.ge.mc) then
                  shat = alfax/ataumi*shat + 9*(1 - alfax/ataumi)/20
     .                 + alfax/pi1*(4*dlog(x/mtau)/15 
     .                 + 21*dlog(alfax/ataumi)/320
     .                 + (11 - 24*zeta3)/900*(atau2 - ax2))
               else
                  call alfahat(mc,dgamma,alfamc)
                  shat = alfamc/ataumi*shat + 9*(1 - alfamc/ataumi)/20
     .                 + alfamc/pi1*(4*dlog(mc/mtau)/15 
     .                 + 21*dlog(alfamc/ataumi)/320
     .                 + (11 - 24*zeta3)/900*(atau2 - ac2))
                  call alfahat(mc - 1.d-14,dgamma,acmins)
                  shat = acmins/alfamc*shat + 3*(1 - acmins/alfamc)/8
                  if (x.ge.mcbar) then
                     shat = alfax/acmins*shat + (1 - alfax/acmins)/2 
     .                    + alfax/pi1*(dlog(x/mc)/3
     .                    + 5*dlog(alfax/acmins)/48)
                  else
                  call alfahat(mcbar,dgamma,amcbar)
                     shat = amcbar/acmins*shat + (1 - amcbar/acmins)/2 
     .                    + amcbar/pi1*(dlog(mcbar/mc)/3 
     .                    + 5*dlog(amcbar/acmins)/48)
                     if (x.ge.msbar) then
                        shat = shat + dlog(mcbar/x)*((1/2.d0 - shat)
     .                       * (dah3mc - dahad2)/dlog(mcbar/msbar) 
     .                       + alfax/3/pi1*(1 - 4*shat)
     .                       * (1 + 3*alfax/4/pi1))
                     else
                        call alfahat(msbar,dgamma,amsbar)
                        shat = shat + (1/2.d0 - shat)*(dah3mc - dahad2) 
     .                       + amsbar/3/pi1*(1 - 4*shat)
     .                       * dlog(mcbar/msbar)*(1 + 3*amsbar/4/pi1)
                        if (x.ge.mudbar) then
                           shat = shat + dlog(msbar/x)*((9/20.d0 - shat)
     .                          * dahad2/dlog(msbar/mudbar) 
     .                          + alfax/3/pi1*(1 - 4*shat)
     .                          * (1 + 3*alfax/4/pi1))
                        else
                           call alfahat(mudbar,dgamma,audbar)
                           shat = shat + (9/20.d0 - shat)*dahad2
     .                          + audbar/3/pi1*(1 - 4*shat)
     .                          *dlog(msbar/mudbar)*(1 + 3*audbar/4/pi1)
                        if (x.ge.mmu) then
                           shat = 1/4.d0 + (shat - 1/4.d0)*alfax/audbar
                           else
                              call alfahat(mmu,dgamma,alfamu)
                              shat = (1 + (4*shat - 1)*alfamu/audbar)/4
                              call alfahat(mmu - 1.d-14,dgamma,amumin)
                              shat = (1 + (4*shat - 1)*amumin/alfamu)/4
                              if (x.ge.me) then
                                 shat = (1 + (4*shat -1)*alfax/amumin)/4
                              else
                                 call alfahat(me,dgamma,alfame)
                                 shat = (1 + (4*shat-1)*alfame/amumin)/4
                                 call alfahat(me - 1.d-14,dgamma,aemins)
                                 shat = (1 + (4*shat-1)*aemins/alfame)/4
                              endif
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
c         if ((flagmr.eqv..false.).or.(mt.le.mz))
c     .        shat = shat - alpha/pi1/6*(1 - 8*sinhat/3)*(dlog(rattz2)
c     .        *(1 + at + alh*xa2mt0/3) - 13/12.d0*at - 5*alh*xa2mt0/4)
      endif
      
      return
      end
