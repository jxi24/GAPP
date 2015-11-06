      subroutine alfahat(x,dgamma,alfaq2)

C   alfahat is computed at scale x. 

      implicit none
      integer nq
      double precision alphas,as,al,alh,alrun,alfaq2,dgamma,x
      double precision cmsbar,const1,const2,const3,dahad2,dah3mc
      double precision mudbar,msbar,mcbar,cut,cut2,msrun
      double precision mcrun,mbrun,match0,match1,match2,me2,mmu2,mtau2
      double precision c0log,c1log,c2log,c3log,mc2,mb2,logmc2,logmb2
      double precision c0cons,c1cons,c2cons,c3cons,mcbcon,cutmc2,cutmb2
      double precision lecut,lmucut,xa2mt0,xalas2,xasmt0
      include 'common.f'		      

      mudbar = 0.176d0
      msbar  = 0.305d0
      mcbar  = 1.176d0

      dah3mc = dahad3*1.175d0
      dahad2 = dahad3*0.298d0

      xasmt0 = 1.d0
      xalas2 = 1.d0
      xa2mt0 = 1.d0
      if (fasmt0.eqv..false.) xasmt0 = 0.d0
      if (falas2.eqv..false.) xalas2 = 0.d0
      if (fa2mt0.eqv..false.) xa2mt0 = 0.d0

      const1 = 5/3.d0
      const2 = 55/12.d0 - 4*zeta3
      const3 = 34525/864.d0 - 9*zeta2/4 - 715*zeta3/18 + 25*zeta5/3
      cmsbar = 1/6.d0

      match0     = 15/ 4.d0*xa2mt0
      match1     = 13/12.d0
      
      al     = alpha/pi1
      alh    = alfaq2/pi1
      cut    = 1.800d0

      me2    = me**2
      mmu2   = mmu**2
      mtau2  = mtau**2
      cut2   = cut**2
      mc2    = mcrun(cut)**2
      mb2    = mbrun(cut)**2
      cutmc2 = cut2/mc2
      cutmb2 = cut2/mb2
      logmc2 = dlog(cutmc2)
      logmb2 = dlog(cutmb2)

      c0log  =      1/12.d0
      c1log  =      2/135.d0
      c2log  = -    1/5040.d0
      c3log  =      1/127575.d0
      c0cons = -   11/72.d0 - const2/6
      c1cons = -    2/25.d0
      c2cons =   1513/2116800.d0
      c3cons = - 1853/80372250.d0

      mcbcon = logmc2   *(c0cons + c0log*logmc2)
     .       + cutmc2   *(c1cons + c1log*logmc2)
     .       + cutmc2**2*(c2cons + c2log*logmc2)
     .       + cutmc2**3*(c3cons + c3log*logmc2)
     .       + cutmb2   *(c1cons + c1log*logmb2)
     .       + cutmb2**2*(c2cons + c2log*logmb2)
     .       + cutmb2**3*(c3cons + c3log*logmb2)

      if (x.lt.mudbar) then
         alh = al
         if (x.ge.me) then
            alh = al/(1 - al**2*match0/4)
             as = 0.d0
            if (x.lt.mmu) then
               alh = alrun(1,0,me,x,alh,as)
            else
               alh = alrun(1,0,me,mmu,alh,as)
               alh = alh/(1 - alh**2*match0/4)
               alh = alrun(2,0,mmu,x,alh,as)
            endif
         endif
      else
             as = alphas(cut)/pi1*xasmt0
          lecut = dlog(cut2/me2)
         lmucut = dlog(cut2/mmu2)
         dgamma = dahad3/al + al*match0/2 + 2/3.d0*(const1 
     .          + (as + al*xa2mt0/4)*(const2 + 2*msrun(cut)**2/cut2)
     .          + as**2*xalas2*(const3 + mcbcon))
     .          + (1/3.d0 + al*xa2mt0/4)*(lecut + lmucut)
     .          + al**2*xa2mt0/24*(lecut**2 + 3*lmucut**2)
         alh = al/(1 - al*dgamma)
         alh = alrun(5,4,cut,mc,alh,as)
          as = alphas(mc)/pi1*xasmt0
         if (x.lt.mc) then
            if (x.ge.mcbar) then
               alh = alrun(5,3,mc,x,alh,as)
            else
               alh = alrun(5,3,mc,mcbar,alh,as)
               if (x.ge.msbar) then
                  alh = alh/(1 + dlog(mcbar/x)*(alh*(4/3.d0 + alh)
     .                + (dah3mc - dahad2)/dlog(mcbar/msbar)))
               else
                  alh = alh/(1 + dlog(mcbar/msbar)*alh*(4/3.d0 + alh)
     .                + (dah3mc - dahad2))
                  alh = alh/(1 + dlog(msbar/x)*(alh*(4/3.d0 + alh)
     .                + dahad2/dlog(msbar/mudbar)))
               endif
            endif
         else
            alh = alh/(1 - 4*alh/9*(as*(match1 + as*match2(4)) 
     .           + alh*match0/3))
            if (x.lt.mtau) then
               alh = alrun(6,4,mc,x,alh,as)
            else
               alh = alrun(6,4,mc,mtau,alh,as)
               as = alphas(mtau)/pi1*xasmt0
               alh = alh/(1 - alh**2*match0/4)
               if (x.lt.mb) then
                  alh = alrun(7,4,mtau,x,alh,as)
               else
                  alh = alrun(7,4,mtau,mb,alh,as)
                  as = alphas(mb)/pi1*xasmt0
                  alh = alh/(1 - alh/9*(as*(match1 + as*match2(5)) 
     .                + alh*match0/12))
                  if (x.lt.mw) then
                     alh = alrun(8,5,mb,x,alh,as)
                  else
                     alh = alrun(8,5,mb,mw,alh,as)
                     as = alphas(mw)/pi1*xasmt0
                     alh = alh/(1 + alh*cmsbar)
                     if (x.lt.mt.and.flagmr.eqv..true.) then
                        alh = alrun(9,5,mw,x,alh,as)
                     else
                        alh = alrun(9,5,mw,mt,alh,as)
                        as = alphas(mt)/pi1*xasmt0
                        alh = alh/(1 - 4*alh/9*(as*(match1 
     .                      + as*match2(6)) + alh*match0/3))
                        alh = alrun(10,6,mt,x,alh,as)
                     endif
                  endif
               endif
            endif
         endif
      endif

      alfaq2 = pi1*alh
      dgamma = 1.d0/al - 1.d0/alh

      return
      end


      double precision function alrun(nf,nq,mu1,mu2,alrun0,as0)

      implicit none
      integer nf,nq
      double precision mu1,mu2,alrun0,as0,l,abetal,ll,as0eff
      double precision beta0(10),beta1(10),delta1(10),singlt(10)
      double precision fbeta0,fbeta1,fbeta2,xalas2,x4lqcd
      double precision b0_1,b1_1,b0_3,c1,c2,del1,del2,del3
      include 'common.f'      
  
      data beta0 /  3,  6, 10, 11, 12, 16, 19, 20, 20, 24/
      data beta1 / 27, 54, 70, 71, 72, 88,115,116,116,132/
      data delta1/  0,  0,  4,  5,  6, 10, 10, 11, 11, 15/
      data singlt/  0,  0,  2,  1,  0,  2,  2,  1,  1,  3/

           l = 2*dlog(mu2/mu1)
        b0_1 = - beta0(nf)/9.d0
        b1_1 = - beta1(nf)/108.d0
        del1 = - delta1(nf)/9.d0
        del2 = del1*(125/48.d0 - 11/72.d0*nq)
        del3 = del1*(10487/1728.d0 + 55*zeta3/18 -
     .              (707/864.d0 + 55*zeta3/54)*nq - 77/3888.d0*nq**2)
     .       - 10*singlt(nf)**2/27.d0*(11/144.d0 - zeta3/6)
        b0_3 = fbeta0(nq)
        c1   = fbeta1(nq)/b0_3
        c2   = fbeta2(nq)/b0_3
      abetal = as0*b0_3*l
          ll = dlog(1 + abetal)
      as0eff = as0/(1 + abetal)

      if (nf.ge.9) b0_1 = b0_1 + 7/4.d0

      xalas2 = 1.d0
      if (falas2.eqv..false.) xalas2 = 0.d0
      x4lqcd = 1.d0
      if (f4lqcd.eqv..false.) x4lqcd = 0.d0
      if (fasmt0.eqv..false.)   del1 = 0.d0
      if (fa2mt0.eqv..false.)   b1_1 = 0.d0

      alrun = 1/(1/alrun0 + b0_1*l + b1_1/b0_1*dlog(1 + alrun0*b0_1*l)
     .      + del1/b0_3*ll + xalas2*as0eff/b0_3*(del2*abetal 
     .      + del1*c1*(ll - abetal)) + x4lqcd*as0eff**2/b0_3
     .      * ((del3 - del2*c1)*abetal + del2*c1*ll - del1*(c1*ll)**2/2
     .      + (del3 - del2*c1 + del1*(c1**2 - c2))*abetal**2/2))

c  For estimate of O(alpha alpha_s^4) contribution:

c      alrun = 1/(1/alrun0 + b0_1*l + b1_1/b0_1*dlog(1 + alrun0*b0_1*l)
c     .      + del1/b0_3*ll + xalas2*as0eff/b0_3*(del2*abetal 
c     .      + del1*c1*(ll - abetal)) + x4lqcd*as0eff**2/b0_3
c     .      * ((del3 - del2*c1)*abetal + del2*c1*ll - del1*(c1*ll)**2/2
c     .      + (del3 - del2*c1 + del1*(c1**2 - c2))*abetal**2/2)
c     .      + as0eff**3/b0_3*27.d0*(abetal + abetal**2 + abetal**3/3)) 

      return
      end


      double precision function match2(nq)

      implicit none
      integer nq
      double precision normq2(6)
      include 'common.f'

C   normq2 is the sum over the Q^2 of the (nq - 1) lighter (massless) quarks,
C   normalized to the Q^2 of the (de)coupling quark.

      data normq2/ 0.d0, 4.d0, 5.d0, 1.5d0, 10.d0, 2.75d0/

C  K.G. Chetyrkin, J.H. K\"uhn, and M. Steinhauser, \npb{482}{213}(1996).
C  2nd line from the Ph. D. thesis of M. Steinhauser for the
C  double bubble diagram with inner massive quark (the one to be matched).

      match2 = 655*zeta3/144 - 3847/864.d0 + nq*361/1296.d0
     .       + normq2(nq)*295/1296.d0

      if (falas2.eqv..false.) match2 = 0.d0

      return
      end
