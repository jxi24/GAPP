      subroutine polemasses(nf,mpole)

C   Running are masses are QCD running masses, but QED pole masses.
C   To convert to QED running masses, a factor 
C   1 - alpha/pi (1 + 3/4 ln mu^2/m^2) needs to be included. 

      implicit none
      double precision msrun,mcrun,mbrun,mtrun,ln2,ln22,ln24
      double precision alphas,as,ratbt,ratct,ratst,ratcb,ratsb,ratsc
      double precision mpole,delta,x,ablm,bblm,cblm,mublm
      integer nf
      include 'common.f'

C   compute quark pole masses in terms of the scale invariant masses 
C   mq := mq(mq); however, these formulae are already in the non-
C   pertubative domain (for the top they are OK) and they do not
C   agree with the one given by K"uhn.

      delta(x) = x*(3/4.d0*zeta2 - 0.597d0*x + 0.230d0*x**2)

      ln2  = dlog(2.d0)
      ln22 =  ln2*ln2
      ln24 = ln22*ln22

      if (nf.lt.4.or.nf.gt.6) then
         print*, 
     .   'Note: pole mass is not defined for this quark flavor!'
         print*,'STOP.'
         stop
      endif
      
      if (nf.eq.4) then
         as = alphas(mc)/pi1
         ratsc = msrun(mc)/mc
         mpole = mc*(1 + 4*as/3 + as**2/3*(779/32.d0 - zeta3/2
     .         + zeta2*(3 + 2*ln2) + 4*delta(ratsc)) 
     .         + as**3/27*(5784469/3456.d0 - (1453*zeta3 - 1975*zeta5)/8
     .         + 488501*zeta2/240 + 185*zeta4/16 - 4317*zeta2*zeta3/8
     .         - (641 + 32*ln2)*ln2*zeta2 - 49*ln24/6 - 196*pol4hf))
      else if (nf.eq.5) then
         as = alphas(mb)/pi1
         ratcb = mcrun(mb)/mb
         ratsb = msrun(mb)/mb
         mpole = mb*(1 + 4*as/3 + as**2/3*(2195/96.d0 - zeta3/2
     .         + zeta2*(2 + 2*ln2) + 4*(delta(ratcb) + delta(ratsb))) 
     .         + as**3/27*(4922965/3456.d0 - 495*zeta3/2 + 1975*zeta5/8
     .         + 439961*zeta2/240 + 1405*zeta4/16 - 4317*zeta2*zeta3/8 
     .         - (663 + 28*ln2)*ln2*zeta2 - 47*ln24/6 - 188*pol4hf))
      else
         ratbt = mbrun(mt)/mt
         ratct = mcrun(mt)/mt
         ratst = msrun(mt)/mt
         if (flgblm.eq.0) then
            as = alphas(mt)/pi1
         mpole = (1 + 4*as/3 + as**2/6*(2053/48.d0 + zeta2*(2 + 4*ln2) 
     .         - zeta3 + 8*(delta(ratbt) + delta(ratct) + delta(ratst)))
     .         + as**3/27*(453365/384.d0 - 2451*zeta3/8 + 1975*zeta5/8
     .         + 394541*zeta2/240 + 2625*zeta4/16 - 4317*zeta2*zeta3/8 
     .         - (685 + 24*ln2)*ln2*zeta2 - 15*ln24/2 - 180*pol4hf))*mt
         else if (flgblm.eq.1) then
            as = alphas(mt)/pi1
         mpole = 1/(1 - 4*as/3 - as**2/6*(1541/48.d0 + zeta2*(2 + 4*ln2) 
     .         - zeta3 + 8*(delta(ratbt) + delta(ratct) + delta(ratst)))
     .         - as**3/27*(280853/384.d0 - 2355*zeta3/8 + 1975*zeta5/8
     .         + 388781*zeta2/240 + 2625*zeta4/16 - 4317*zeta2*zeta3/8
     .         - (733 + 24*ln2)*ln2*zeta2 - 15*ln24/2 - 180*pol4hf))*mt
         else
            ablm =    4/3.d0
            bblm = 307/32.d0 + pi2/3 + 2*dlog(2.d0)*zeta2/3 - zeta3/6
     .           + 4/3.d0*(delta(ratbt) + delta(ratct) + delta(ratst))
            cblm = 71/144.d0 + zeta2/3         
            mublm = dexp( - 3*cblm/ablm)*mt
            as = alphas(mublm)/pi1
            mpole = mt*(1 + ablm*as + (bblm - 33/2.d0*cblm)*as**2)
         endif
      endif

      return
      end


      double precision function msrun(x)

      implicit none
      double precision x,a1,ax,alphas,mrun,amatch,mmatch
      include 'common.f'

      a1 = alphas(1.d0)/pi1
      amatch = 11/72.d0
      mmatch = 89/432.d0

      if (x.ge.mc) then
         ax = alphas(mc)/pi1
         msrun = mrun(3,a1,ax*(1 + amatch*ax**2),ms)
         msrun = msrun/(1 + mmatch*ax**2)
         a1 = ax
         if (x.ge.mb) then
            ax = alphas(mb)/pi1
            msrun = mrun(4,a1,ax*(1 + amatch*ax**2),msrun)
            msrun = msrun/(1 + mmatch*ax**2)
            a1 = ax
            if (x.ge.mt) then
               ax = alphas(mt)/pi1
               msrun = mrun(5,a1,ax*(1 + amatch*ax**2),msrun)
               msrun = msrun/(1 + mmatch*ax**2)
               a1 = ax
               ax = alphas(x)/pi1                  
               msrun = mrun(6,a1,ax,msrun)
            else
               ax = alphas(x)/pi1                  
               msrun = mrun(5,a1,ax,msrun)
            endif
         else
            ax = alphas(x)/pi1                  
            msrun = mrun(4,a1,ax,msrun)
         endif            
      else
         ax = alphas(x)/pi1 
         msrun = mrun(3,a1,ax,ms)
      endif
      
      return
      end


      double precision function mcrun(x)

      implicit none
      double precision x,a1,ax,alphas,mrun,amatch,mmatch
      include 'common.f'

      a1 = alphas(mc)/pi1
      amatch = 11/72.d0
      mmatch = 89/432.d0

      if (x.ge.mb) then
         ax = alphas(mb)/pi1
         mcrun = mrun(4,a1,ax*(1 + amatch*ax**2),mc)
         mcrun = mcrun/(1 + mmatch*ax**2)
         a1 = ax
         if (x.ge.mt) then
            ax = alphas(mt)/pi1                  
            mcrun = mrun(5,a1,ax*(1 + amatch*ax**2),mcrun)
            mcrun = mcrun/(1 + mmatch*ax**2)
            a1 = ax
            ax = alphas(x)/pi1                  
            mcrun = mrun(6,a1,ax,mcrun)
         else
            ax = alphas(x)/pi1                  
            mcrun = mrun(5,a1,ax,mcrun)
         endif
      else
         ax = alphas(x)/pi1 
         if (x.lt.mc) then
            mcrun = mc*(1 + mmatch*a1**2)
            mcrun = mrun(3,a1*(1 + amatch*a1**2),ax,mcrun)
         else
            mcrun = mrun(4,a1,ax,mc)
         endif
      endif
      
      return
      end


      double precision function mbrun(x)

      implicit none
      double precision x,a1,ax,alphas,mrun,amatch,mmatch
      include 'common.f'

      a1 = alphas(mb)/pi1
      amatch = 11/72.d0
      mmatch = 89/432.d0

      if (x.ge.mb) then         
         if (x.ge.mt) then
            ax = alphas(mt)/pi1
            mbrun = mrun(5,a1,ax*(1 + amatch*ax**2),mb)
            mbrun = mbrun/(1 + mmatch*ax**2)
            a1 = ax
            ax = alphas(x)/pi1                  
            mbrun = mrun(6,a1,ax,mbrun)
         else
            ax = alphas(x)/pi1                  
            mbrun = mrun(5,a1,ax,mb)
         endif
      else
         mbrun = mb*(1 + mmatch*a1**2)
         if (x.lt.mc) then
            ax = alphas(mc)/pi1
            mbrun = mrun(4,a1*(1 + amatch*a1**2),ax,mbrun)
            a1 = ax
            ax = alphas(x)/pi1                  
            mbrun = mbrun*(1 + mmatch*a1**2)
            mbrun = mrun(3,a1*(1 + amatch*a1**2),ax,mbrun)
         else
            ax = alphas(x)/pi1 
            mbrun = mrun(4,a1*(1 + amatch*a1**2),ax,mbrun)
         endif
      endif

      return
      end


      double precision function mtrun(x)

      implicit none
      double precision x,a1,ax,alphas,mrun,amatch,mmatch
      include 'common.f'

      a1 = alphas(mt)/pi1
      amatch = 11/72.d0
      mmatch = 89/432.d0

      if (x.lt.mt) then
         mtrun = mt*(1 + mmatch*a1**2)         
         if (x.lt.mb) then
            ax = alphas(mb)/pi1
            mtrun = mrun(5,a1*(1 + amatch*a1**2),ax,mtrun)
            a1 = ax
            mtrun = mtrun*(1 + mmatch*a1**2)
            if (x.lt.mc) then
               ax = alphas(mc)/pi1
               mtrun = mrun(4,a1*(1 + amatch*a1**2),ax,mtrun)
               a1 = ax
               mtrun = mtrun*(1 + mmatch*a1**2)
               ax = alphas(x)/pi1                  
               mtrun = mrun(3,a1*(1 + amatch*a1**2),ax,mtrun)
            else
               ax = alphas(x)/pi1                  
               mtrun = mrun(4,a1*(1 + amatch*a1**2),ax,mtrun)
            endif
         else
            ax = alphas(x)/pi1                  
            mtrun = mrun(5,a1*(1 + amatch*a1**2),ax,mtrun)
         endif
      else
         ax = alphas(x)/pi1                  
         mtrun = mrun(6,a1,ax,mt)
      endif

      return
      end


      double precision function mrun(nf,a1,ax,mrun0)

      double precision a1,ax,mrun0,gamma1,gamma2,fac1,fac2
      double precision beta0,beta1,beta2,fbeta0,fbeta1,fbeta2
      integer nf
      include 'common.f'      
      
      if (a1.eq.ax) then
         mrun = mrun0
         return
      endif

      beta0 = fbeta0(nf)
      beta1 = fbeta1(nf)
      beta2 = fbeta2(nf)
      
      gamma1 = (202/3.d0 - 20/9.d0*nf)/16
      gamma2 = (1249 - (2216/27.d0 + 160/3.d0*zeta3)*nf 
     .               -   140/81.d0*nf**2)/64

      fac1 = gamma1/beta0 - beta1/beta0**2
      fac2 = gamma2/beta0 - beta2/beta0**2
      
      mrun = mrun0*(ax/a1)**(1.d0/beta0)*exp(fac1*(ax - a1) + 
     .            (fac2 - beta1/beta0*fac1)*(ax**2 - a1**2)/2)

      return
      end
