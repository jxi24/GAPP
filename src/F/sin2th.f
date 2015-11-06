      subroutine sin2thetaw

C   following G. Degrassi,   S. Fanchiotti and A. Sirlin, NPB 351 (1991)  49,
C         and S. Fanchiotti, B. Kniehl     and A. Sirlin, PRD  48 (1993) 307.
C  
C   We improve on these references by using MS-bar fermion mass definitions.
C
C   For delta r_W we take m_c(M_W) rather than m_c(M_Z), since we 
C   checked that then no ln(cos2tw) terms arise in O(alpha alpha_s m_c^2). 
C   We take m_b(m_t), since this is the scale for the t-b doublet, and we have
C   checked that a high scale m_b eliminates all fermion mass singularities in 
C   O(alpha alpha_s m_b^2). We then neglect the very complicated expression 
C   in that order. Accordingly, alpha_s(M_W) is chosen for the light quarks,
C   and alpha_s(m_t) for the t-b doublet. Here, we also neglect the small 
C   O(alpha^2 m_t^2) corrections.
C
C   For delta r_Z we take m_c(M_Z), m_b(M_Z) and alpha_s(M_Z) except for the
C   top contribution, for which we take alpha_s(m_t).
C
C   We did not yet include O(alpha alpha_s^2) corrections subleading
C   in m_t^2. They are computed for delta_r and delta_kappa (on-shell) by
C   K.G. Chetyrkin, J.H. K"uhn and M. Steinhauser in SLAC-PUB-95-6851 
C   (submitted to PRL).

      implicit none
      integer i,f
      double precision alphas,az,aw,at,al,alh,dgamma,const,cqcd1,cqcd2
      double precision mbrun,mcrun,t1,t2,t3,c(2,3),xa2mt0,xasmt0,rootc
      double precision zt,b1,b2,b3,b4,b5,b6,b0zww,b0wwz,b0zhz,b0whw
      double precision rtauz2,ratcz2,ratbz2,rtauw2,ratcw2,ratbw2
      double precision i1dfs,i2dfs,hms,hmsw,fftj0,fftj1,gftj,f1bk,b0
      double precision a0,delrho,xt,fftj00,fftj10,gftj0,weight,ddilog
      double precision lambda,lambd1,contrf,contmb,conqcd,dfks,dfksh
      double precision sin2tw,cos2tw,cos4tw,sin4th,cos4th,c6,s4c2,cfunc
      complex*16 dkapse
      include 'common.f'

C   Relevant Coefficients when using on-shell top-mass definition;
C   in that case also the constant at and at*mt2/mz2 coefficients 
C   (for the latter cqcd1 ---> - 1/2 - zeta2) must be changed and also cqcd2.
C
C     c(1,1) =      82/     81.d0
C     c(1,2) =     449/   2700.d0
C     c(1,3) =   62479/1984500.d0
C     c(2,1) =     689/   1620.d0
C     c(2,2) =    1691/  37800.d0
C     c(2,3) =   49213/7938000.d0

      c(1,1) =     194/    405.d0
      c(1,2) =     983/  18900.d0
      c(1,3) =   12079/1984500.d0
      c(2,1) =     257/   1620.d0
      c(2,2) =     251/  37800.d0
      c(2,3) = -  1187/7938000.d0

          al = alpha/pi1
          az = alphas(mz)/pi1
          at = alphas(mt)/pi1
          xt = gf*mt2/2/pi2/dsqrt(2.d0)
       const = (55/12.d0 - 4*zeta3)/3
       cqcd1 = zeta2 - 3/2.d0
       cqcd2 = 2.9772d0
      rtauz2 = (mtau     /mz)**2
      ratcz2 = (mcrun(mz)/mz)**2
      ratbz2 = (mbrun(mz)/mz)**2

      xa2mt0 = 1.d0
      xasmt0 = 1.d0
      if (fa2mt0.eqv..false.) xa2mt0 = 0.d0
      if (fasmt0.eqv..false.) xasmt0 = 0.d0

C   first order iteration:

      a0  = dsqrt(pi1*alpha/dsqrt(2.d0)/gf)
      mw2 = mz2/2*(1 + dsqrt(1 - 4*a0**2/mz2))

      if (flagzp.ne.0) mw2 = mz02/2*(1 + dsqrt(1 - 4*a0**2/mz02))
      mw  = dsqrt(mw2)

      alphat = alpha
      call alfahat(mz,dgamma,alphat)

      delrwh = al*dgamma
      rhohat = 1 + 3*xt/4
      mw2 = mz2*rhohat/2*(1 + dsqrt(1 - 4*a0**2/mz2/rhohat/(1- delrwh)))
      if (flagzp.ne.0) mw2 = mz02*rhohat/2
     .     *(1 + dsqrt(1 - 4*a0**2/mz02/rhohat/(1 - delrwh)))
      mw  = dsqrt(mw2)
      cos2tw = mw2/mz2
      sin2tw = 1 - cos2tw
      coshat = cos2tw/rhohat
      sinhat = 1 - coshat

      do 10 f = 0, 9
         v(f) = (i3(f) - 2*q(f)*sinhat)/2/dsqrt(sinhat*coshat)
         a(f) =  i3(f)/2/dsqrt(sinhat*coshat)
 10   continue

C   loop iteration:

      do 100 i = 1,3

      call alfahat(mz,dgamma,alphat)

         alh = alphat/pi1
          aw = alphas(mw)/pi1
      rtauw2 = (mtau     /mw)**2
      ratcw2 = (mcrun(mw)/mw)**2
      ratbw2 = (mbrun(mt)/mw)**2

      sin4th = sinhat**2
      cos4tw = 1/cos2tw**2
        s4c2 = sin4th/coshat
       cfunc = (1/cos2tw - cos4tw/6)/4
       rootc = dsqrt(4*cos2tw - 1)
        hmsw = hms(mh2/mw2)

C   Vertex + box corrections and decoupling effects in delrwh:

      delrwh = al*dgamma + al/4/sinhat*(6 + dlog(cos2tw)/sin2tw*
     .     (7/2.d0 - 5*sin2tw/2 - sinhat*(5 - 3*rhohat/2)))

C  functions (A.10) and (A.7) in Degrassi, Fanchiotti and Sirlin.

      i1dfs = (21*coshat/2 + s4c2/2 - sinhat)/18 - (lambda(1/cos2tw)- 1)
     .     *(cos2tw*(10*coshat - 10*s4c2/3 - 4*sinhat/3) 
     .     + 13*coshat/2 + sinhat/3 - s4c2/6)

      i2dfs = rootc/cos2tw*datan(rootc)*(s4c2 - 1/cos2tw/3 + cos4tw/12 - 
     .     coshat*(5 + 11/cos2tw/3 - 2*cos4tw/3)) + dlog(cos2tw)/cos2tw*
     .     (coshat*(5/cos2tw/2 - cos4tw/3 + 7/2.d0) + cfunc - s4c2/2) + 
     .     coshat*(101/9.d0 + 6/cos2tw - 2*cos4tw/3) + 
     .     1/9.d0 + 2*(cfunc - s4c2)

C   Bosonic contribution to delrwh:

      delrwh = delrwh + al/4/sinhat*(hmsw + i2dfs - 3/4.d0*
     .     (rathz2*dlog(rathz2) - cos2tw*dlog(cos2tw))/(rathz2 - cos2tw)
     .     + dlog(cos2tw)*(-(4*coshat + 1/4.d0)/cos2tw - 20/3.d0 +
     .     4*sinhat + s4c2 - cos2tw/sin2tw*(4*coshat + 1/4.d0 - s4c2))
     .     - 7/4.d0 + 92*sinhat/9 - (2*coshat + 1/8.d0)/cos2tw + s4c2)

C   Fermionic contribution to delrwh:

      rattw2 = mt2/mw2
      contrf = 1/4.d0/sinhat*(4*dlog(cos2tw) - 20/3.d0 +
     .     dlog(rattw2) + rattw2/4*(1 + 2*rattw2) + 
     .     (rattw2 - 1)**2*(1 + rattw2/2)*dlog(dabs(1 - 1/rattw2)))
      delrwh = delrwh + al*contrf

C   O(alpha alpha_s) corrections to delrwh; top quark term according to
C   B. A. Kniehl, NPB 347 (1990) 86: 

      if (fasmt0.eqv..true.) then
         conqcd = aw/2/sinhat*(dlog(cos2tw) - 3*const) +
     .        at/4/sinhat*(dlog(rattz2) - 3*const -
     .        4*rattw2*(f1bk(mw2/mt2) - 23/16.d0 + zeta2/2 + 3*zeta3/2))

C   Extra term from conversion to MS-bar top mass:

         conqcd = conqcd + at/2/sinhat*rattw2*
     .        (1 + 2*rattw2 - 2*(1 - rattw2**2)*dlog(dabs(1 - mw2/mt2)))

         delrwh = delrwh + al*conqcd
      endif

C   Bosonic contribution to rhohat:

      delrho = alh/4/sinhat/coshat*(cos2tw*hmsw - hms(rathz2) +
     .     cos2tw*i2dfs - coshat*i1dfs + coshat**2/2 + sin4th/6 -
     .     cos2tw*(1/2.d0 - 83/9.d0*sinhat) - dlog(cos2tw)*(1 +
     .     cos2tw*(7/6.d0 + 3*sinhat - s4c2) - 3*sin4th + 7*sinhat/3))

C   Fermionic contribution to rhohat:

      lambd1 = lambda(1/rattz2) - 1
      delrho = delrho + alh*rhohat*contrf +
     .     3*alh/4/sinhat/coshat*rattz2*(1/4.d0 + lambd1) +
     .     alh/2/sinhat/coshat*((7/4.d0 - 10*sinhat/3 + 40*sin4th/9)*
     .     5/3.d0 - (1 + (1 - 8*sinhat/3)**2)/8*
     .     (dlog(rattz2) + 1/3.d0 + 2*(1 + 2*rattz2)*lambd1))

C   Leading O(alpha alpha_s mt**2) corrections to rhohat:

      if ((fasmt2.eqv..true.).and.(ffermi.eqv..false.)) delrho = delrho 
     .     - alh*at/4/sinhat/coshat*cqcd1*rattz2

C   Remaining O(alpha alpha_s) corrections to rhohat:

      if (fasmt0.eqv..true.) then
         t1 = 1.d0/rattz2
         t2 = 1.d0/rattz2**2
         t3 = 1.d0/rattz2**3
         delrho = delrho + alh*rhohat*conqcd + alh*az/4/sinhat/coshat*
     .        3*const*(5/2.d0 - 14*sinhat/3 + 44*sin4th/9) +
     .        alh*at/4/sinhat/coshat*( - dlog(rattz2)*
     .        (1/2.d0 - 4*sinhat/3 + 16*sin4th/9) +
     .        (1/4.d0 - 4*sinhat/3 + 16*sin4th/9)*
     .        (13/12.d0 + c(1,1)*t1 + c(1,2)*t2 + c(1,3)*t3) + 1/4.d0*
     .       (-29/36.d0 + c(2,1)*t1 + c(2,2)*t2 + c(2,3)*t3))
      endif

C   Leading O(alpha alpha_s**2 mt**2) correction to rhohat:
      
      if ((fas2mt.eqv..true.).and.(ffermi.eqv..false.)) delrho = delrho 
     .     - alh*at**2/4/sinhat/coshat*rattz2*cqcd2
      
C   Leading O(alpha**2 mt**4) term with full M_H dependence as comupted by
C   R. Barbieri, M. Beccaria, P. Ciafaloni, G. Curci and A. Vicere, 
C   PLB 288 (1992) 95 and NPB 409 (1993) 105; we use the expressions from
C   J. Fleischer, O.V. Tarasov and F. Jegerlehner, PLB 319 (1993) 249.

      rho2 = 0.d0
      if ((fa2mt4.eqv..true.).and.(fa2mt2.eqv..false.)) then
         b1 = 1/ratth2
         b2 = b1**2
         rho2 = 25 - 4*b1 + (b2 - 12*b1 - 12)/2*dlog(b1)
     .        + (b1 - 2)/2/b1*pi2 + (b1 - 4)/2*dsqrt(b1)*gftj(b1)
     .        - (3*b2 - 12*b1 + 15 - 6/b1)*fftj0(b1) 
     .        + 3*(b2 - 6*b1 + 10)*fftj1(b1)
      endif

C   O(alpha**2 mt**4) +  O(alpha**2 mt**2) corrections to delrho;
C   approximation for mh/mt < 3.8 according to e-mail 493 on scipp 
C   by P. Gambino 3/30/96 (uses sinhat = 0.2314)
C
C      rho2 = - 15.642064d0 + 0.036381841d0*mt +
C     .     ( -  9.95272d0  + 0.0180877d0  *mt)*ratht +
C     .     (    5.68703d0  - 0.0156807d0  *mt)*ratht**2 +
C     .     ( -  1.64687d0  + 0.00536879d0 *mt)*ratht**3 +
C     .     (    0.185188d0 - 0.000646014d0*mt)*ratht**4 +
C     .     (    2.30111d0  - 0.013429d0   *mt)*dsqrt(ratht)
C  
C   we use light and heavy Higgs expansions as given in e-mail 497 on scipp 
C   by P. Gambino 3/31/96, which differs in log(zt) terms from the paper
C   G. Degrassi, P. Gambino and A. Vicini, hep-ph/9603374.
C
C   For mh >     mt, we use the heavy Higgs expansion;
C   for mh < 0.6*mt, we use the light Higgs expansion;
C   in between we use a linear interpolation. 

      if (fa2mt2.eqv..true.) then
          if (mh.le.0.6d0*mt) then
              zt = 1/rattz2
              b1 = 1/ratth2
              b2 = b1*b1
          cos4th = coshat**2

           b0zww = b0(mz,mz2,mw2,mw2)
           b0wwz = b0(mz,mw2,mw2,mz2)
           b0zhz = b0(mz,mz2,mh2,mz2)
           b0whw = b0(mz,mw2,mh2,mw2)
                 
            rho2 = 19 - 53*b1/3 + 3*pi1/2*dsqrt(b1)**3 + 8*b2/9/zt
     .           - 5*b2/9/coshat/zt - 4*dsqrt(b1)*pi1
     .           + (845 - 9/coshat + 427*coshat - 366*cos4th)/27*zt
     .           + (54*b1 - 54 - 119*zt + 44*coshat*zt)/27*pi2
     .           + (34 - 116*coshat + 64*cos4th)*zt*4*pi1/27*dsqrt(b1)
     .           +            (32*b1 - 8*b2/zt -   96*zt)/9*b0zhz
     .           +  zt/coshat*(1 + 20*coshat - 24*cos4th)/3*b0wwz
     .           -       2*zt*(1 + 18*coshat - 16*cos4th)/3*b0zww
     .           - 5*(4*b1 - b2/coshat/zt - 12*coshat*zt)/9*b0whw
     .           - dlog(coshat)*
     .             (5*b1/zt + 3 + 32*coshat + 48*cos4th)/9*zt 
     .           + (5*b2/coshat - 8*b2 - 18*b1*zt)/9/zt*dlog(b1)
     .           + ( - b1 + 8*b2/9/zt - 5*b2/9/coshat/zt - zt/9 
     .           + 104*coshat*zt/3 - 104*cos4th*zt/9)*dlog(zt)

         else if (mh.ge.mt) then

              zt = 1/rattz2
              b1 = 1/ratth2
              b2 = b1*b1
              b3 = b1*b2
              b4 = b1*b3
              b5 = b1*b4
              b6 = b1*b5
          cos4th = coshat**2
              c6 = coshat**3

          fftj00 = fftj0(b1) 
          fftj10 = fftj1(b1) 
           gftj0 =  gftj(b1) 

           b0zww = b0(mz,mz2,mw2,mw2)
           b0wwz = b0(mz,mw2,mw2,mz2)
                
            rho2 = 25 - 4*b1 + (b2 - 12*b1 - 12)/2*dlog(b1)
     .           + (1 - 2/b1)/2*pi2 + (b1 - 4)/2*dsqrt(b1)*gftj0
     .           - (3*b2 - 12*b1 + 15 - 6/b1)*fftj00 
     .           + 3*(b2 - 6*b1 + 10)*fftj10 + zt*(( - 1776*cos4th
     .           + (   72 - 6250*coshat - 3056*cos4th + 3696*c6)*b1
     .           + ( - 18 + 1283*coshat + 1371*cos4th - 1436*c6)*b2  
     .           + (          68*coshat -  124*cos4th +  128*c6)*b3)/
     .             54/coshat/(b1 - 4)/b1 + pi2/27/b2*
     .             ( - 37*coshat + 6*coshat*b1 - 119*b2 + 56*coshat*b2)
     .           + ( - 2/3.d0 - 12*coshat   + 32*cos4th/3)*b0zww
     .           + (  20/3.d0 +  1/coshat/3 -  8*coshat  )*b0wwz
     .           + (17 - 58*coshat + 32*cos4th)*(4 - b1)*
     .             dsqrt(b1)*gftj0/27 + 20*sinhat*(4 - b1)*
     .             gftj0/3/dsqrt(b1)/b1 + ddilog(1 - b1)*
     .             (37 - 6*b1 - 12*b2 - 22*b3 + 9*b4)*2*coshat/9/b2
     .           - (1 + 14*coshat + 16*cos4th)*dlog(coshat)/3)
            rho2 = rho2
     .           + zt*((11520 - 15072*coshat 
     .           - ( 7170 -  8928*coshat -  768*cos4th)*b1 
     .           + ( 3411 -  7062*coshat + 3264*cos4th)*b2 
     .           - ( 1259 -  3547*coshat + 2144*cos4th)*b3 
     .           + (  238 -   758*coshat +  448*cos4th)*b4 
     .           - (   17 -    58*coshat +   32*cos4th)*b5)*
     .             dlog(b1)/27/(b1 - 4)**2/b1 
     .           - (   81 -   362*coshat +  104*cos4th)/9*dlog(zt)
     .           - 2*fftj10/9/(b1 - 4)**2/b2*(3840*sinhat
     .           - ( 4310 -  4224*coshat -  256*cos4th)*b1 
     .           + ( 1706 -  1312*coshat -  320*cos4th)*b2 
     .           - (  315 +   476*coshat -   64*cos4th)*b3
     .           + (24 + 454*coshat)*b4 - 112*coshat*b5 + 9*coshat*b6))

            if (b1.le.4) rho2 = rho2 
     .                        - zt*20*pi1/3*sinhat*dsqrt(4/b1 - 1)**3

        else

              zt = 1/rattz2
              b1 = 0.3600d0
              b2 = 0.1296d0
          cos4th = coshat**2 
          weight = (mh/mt - 0.6d0)/0.4d0

          fftj10 = fftj1(1.d0) 
           gftj0 =  gftj(1.d0) 

           b0zww = b0(mz,mz2,mw2,mw2)
           b0wwz = b0(mz,mw2,mw2,mz2)
           b0zhz = b0(mz,mz2,b1*mt2,mz2)
           b0whw = b0(mz,mw2,b1*mt2,mw2)

            rho2 = weight*(21 - pi2/2 - 3/2.d0*gftj0 + 15*fftj10
     .           + ((1633 - 18/coshat + 1195*coshat - 796*cos4th)/54
     .           + pi2/27*(25*coshat - 119) - 20*pi1*sinhat*dsqrt(3.d0)
     .           + ( - 2/3.d0 - 12*coshat   + 32*cos4th/3)*b0zww
     .           + (  20/3.d0 +  1/coshat/3 -  8*coshat  )*b0wwz
     .           - (1 + 14*coshat + 16*cos4th)*dlog(coshat)/3 
     .           - ( 81 - 362*coshat + 104*cos4th)/9*dlog(zt)
     .           + (197 - 238*coshat + 32*cos4th)*gftj0/9
     .           - fftj10*(26*sinhat - 8/3.d0))*zt)
            rho2 = rho2
     .           + (1 - weight)*(19 - 53*b1/3 + 3*pi1/2*dsqrt(b1)**3 
     .           + 8*b2/9/zt - 5*b2/9/coshat/zt - 4*dsqrt(b1)*pi1
     .           + (845 - 9/coshat + 427*coshat - 366*cos4th)/27*zt
     .           + (54*b1 - 54 - 119*zt + 44*coshat*zt)/27*pi2
     .           + (34 - 116*coshat + 64*cos4th)*zt*4*pi1/27*dsqrt(b1)
     .           +            (32*b1 - 8*b2/zt -   96*zt)/9*b0zhz
     .           +  zt/coshat*(1 + 20*coshat - 24*cos4th)/3*b0wwz
     .           -       2*zt*(1 + 18*coshat - 16*cos4th)/3*b0zww
     .           - 5*(4*b1 - b2/coshat/zt - 12*coshat*zt)/9*b0whw
     .           - dlog(coshat)*
     .             (5*b1/zt + 3 + 32*coshat + 48*cos4th)/9*zt 
     .           + (5*b2/coshat - 8*b2 - 18*b1*zt)/9/zt*dlog(b1)
     .           + ( - b1 + 8*b2/9/zt - 5*b2/9/coshat/zt - zt/9 
     .           + 104*coshat*zt/3 - 104*cos4th*zt/9)*dlog(zt))

         endif
      endif
      
      if (ffermi.eqv..false.) delrho = delrho + 
     .     3*(alh/16/sinhat/coshat*rattz2)**2*rho2

C   Replacement in the leading mt**2 and mt**4 terms of 
C   alphat/4/pi/sinhat/coshat by gf*mz**2/2/dsqrt(2)/pi2*rhohat, i.e.
C   alphat/4/pi/sinhat        by gf*mw**2/2/dsqrt(2)/pi2.
C   However, when O(alpha**2 mt**2) corrections are included, this must not
C   be done in the O(alpha mt**2) term.

      if (ffermi.eqv..true.) then
         if (fa2mt2.eqv..false.) then
            delrho = delrho - 3*alh/16/sinhat/coshat*rattz2 
     .             + 3*xt/4*rhohat*(1 + xt*rhohat/4*rho2)
         else
            delrho = delrho + 3*(xt*rhohat)**2/16*rho2
         endif
         if (fasmt2.eqv..true.) delrho = delrho -    at*xt*rhohat*cqcd1
         if (fas2mt.eqv..true.) delrho = delrho - at**2*xt*rhohat*cqcd2
      endif

      if (fla2im.eqv..true.) then
         call deltakappase(dkapse,.true.)
         delrho = delrho + sinhat/coshat*dimag(dkapse)**2
      endif

      if ((flagmr.eqv..true.).and.(mt.gt.mz)) then
           dfks = (1/sinhat/6 - 4/9.d0)*(dlog(rattz2)*(1 + xasmt0*at + 
     .            xa2mt0*alh/3) - xasmt0*13/12*at - xa2mt0*5/4*alh)
          dfksh = dfks/(1 - 8*sinhat/3)
         delrwh = delrwh - al*dfksh
         delrho = delrho*(1 - alh*dfksh) - alh*sinhat/coshat*dfks
      endif

      if (fobliq.eqv..true.) then
         delrwh = delrwh + alpha/4/sinhat*SWpar
         delrho = delrho + alphat
     .          * (Tpar + (cos2tw*SWpar - Spar)/4/sinhat/coshat)
      endif

      delrzh = delrwh - delrho*(1 - delrwh)

C   Subsubleading (only photonic) O(alpha**2) corrections to delrzh:

      if (fa2mt0.eqv..true.) delrzh = delrzh + al*alh/4/sinhat/coshat*
     .    (const*(5/2.d0 - 23*sinhat/3 + 116*sin4th/9) - 17/27.d0 -
     .    (1/2.d0 - 4*sinhat/3 + 16*sin4th/9)*(dlog(rattz2)/3 - 5/4.d0))

C   light fermion mass effects:

      if (flagmf.eqv..true.) then
         contmb = - ((1 + rattw2**2)*dlog(dabs(1 - mw2/mt2)) + rattw2)
         delrwh = delrwh + al/8/sinhat*
     .        (3*ratbw2*(     contmb  - 1/2.d0)
     .       + 3*ratcw2*(dlog(ratcw2) - 1/2.d0)
     .       +   rtauw2*(dlog(rtauw2) - 1/2.d0))

         delrzh = delrzh + al/8/sinhat/coshat*
     .    (3*ratbz2*(dlog(rattz2) - 16*sin4th/9 +  8*sinhat/3 - 1/2.d0)
     .   + 3*ratcz2*(dlog(ratcz2) - 64*sin4th/9 + 16*sinhat/3 - 1/2.d0)
     .   +   rtauz2*(dlog(rtauz2) - 16*sin4th   +  8*sinhat   - 1/2.d0))

         if (fasmt0.eqv..true.) then
            delrwh = delrwh + al*aw/4/sinhat*ratcw2*(31/4.d0 - 
     .           8*zeta2 - 6*zeta3 + dlog(ratcw2) - 3*dlog(ratcw2)**2/2)

            delrzh = delrzh + al*az/4/sinhat/coshat*
     .           (ratcz2*(31/4.d0 - 8*zeta2 - 6*zeta3 +   dlog(ratcz2) - 
     .           3*dlog(ratcz2)**2/2 - 64*sinhat/9*(4*sinhat - 3)) +
     .           ratbz2*(43/4.d0 - 10*zeta2 - 6*zeta3 + 4*dlog(rattz2) - 
     .           3*dlog(rattz2)**2/2 - 32*sinhat/9*(2*sinhat - 3)))
         endif
      endif

      sinhat = (1 - dsqrt(1 - 4*a0**2/mz2/(1 - delrzh)))/2
      if (flagzp.ne.0) sinhat = 
     .     (1 - dsqrt(1 - 4*a0**2/mz02/(1 - delrzh)))/2
      if (sinhat.eq.0.d0) sinhat = 1.d-15
      coshat = 1 - sinhat
         mw2 = a0**2/sinhat/(1- delrwh)
          mw = dsqrt(mw2)
      cos2tw = mw2/mz2
      sin2tw = 1 - cos2tw
      rhohat = cos2tw/coshat

      do 20 f = 0, 9
         v(f) = (i3(f) - 2*q(f)*sinhat)/2/dsqrt(sinhat*coshat)
         a(f) =  i3(f)/2/dsqrt(sinhat*coshat)
 20   continue

 100  continue

C Kai015 ---------------------------------------
C     Corrections to the W mass
C ----------------------------------------------

      if( modtype.eq.1 ) then
      
        mw2 = mw2 + 
     -  (alpha*(fitcph**2 - fits2b)*kkcc*pi1)/
     -  (dsqrt(2.d0)*fitx*gf*(kkcc - kkss)*kkss)

      endif
      
      if( modtype.eq.2 ) then
      
        mw2 = mw2 + 
     -  (alpha*(fitcph**2 - 2*fits2b)*kkcc*pi1)/
     -  (4.*dsqrt(2.d0)*fitx*gf*(kkcc - kkss)*kkss)

      endif

      if( modtype.eq.3 ) then

        mw2 = mw2 + 
     -  (alpha*(-1 + fitcph)**2*pi1)/
     -  (dsqrt(2.d0)*fitx*gf*(kkcc - kkss))

      endif

      if( modtype.eq.4 ) then

        mw2 = mw2 + 
     -  (alpha*(-1 + fitcph)**2*pi1)/
     -  (dsqrt(2.d0)*fitx*gf*(kkcc - kkss))

      endif

      mw = dsqrt(mw2)

C Kai015 ---------------------------------------

      delr = delrwh - coshat/sinhat*(delrwh - delrzh)/
     .     (1 - coshat/sinhat*(delrwh - delrzh)/(1 - delrwh))

      ratzw2 = mz2/mw2
      rathw2 = mh2/mw2
      rattw2 = mt2/mw2

      return 
      end


      double precision function hms(x)
      
      implicit none
      double precision x,xx,a1,a2

      xx = x**2

      if (x.ge.4.d0) then
         a2 = xx/4.d0 - x
         a1 = dsqrt(a2)
         hms = ((1 + a2/2)*(x/2 - a1) + (a1**3 - x**3/8 + 3*x/2)/6)*
     .        dlog(x) + 3*x/8 - xx/12 - 17/9.d0 + 2*a1*(1 + a2/3)*
     .        dlog((a1 + x/2)/(a1 - 1 + x/2))
      else
         hms = x*dlog(x)/4*(3 - x + xx/6) + 3*x/8 - xx/12 - 17/9.d0 +
     .       dsqrt(x*(1 - x/4))*(2 - 2*x/3 + xx/6)*datan(dsqrt(4/x - 1))
      endif

      return
      end


      double precision function f1bk(x)
      
      implicit none
      double precision x,xx,a0,b,ddilog,dtrilog
      include 'common.f'

      xx = x**2

      if (x.le.1.d0) then
         b = dlog(1 - x)
         f1bk = (x - 3/2.d0 + 1/xx/2)*(b/3*(2*ddilog(x) - zeta2) +
     .        b**3/6 - dtrilog(x) - dtrilog( - x/(1 - x))) +
     .        (x + 1/2.d0 - 1/x/2)/3*ddilog(x) +
     .        (x - 3/4.d0 - 3/x/2 + 5/xx/4)/6*b**2 -
     .        (x - 5/2.d0 + 2/x/3 + 5/xx/6)/4*b + zeta3*(x - 3/2.d0) +
     .        zeta2/3*(x - 7/4.d0 - 1/x/2) + 13/12.d0 - 5/x/24
      else
         a0 = dlog(x)
          b = dlog(x - 1)
         f1bk = (x - 3/2.d0 + 1/xx/2)*(dtrilog(1/(1 - x)) +
     .        2*b/3*ddilog(1/(1 - x)) + (zeta2 - b**2/6)*(a0 - b)) +
     .        (x + 1/2.d0 - 1/x/2)/3*(ddilog(1/(1 - x)) - a0*b) +
     .        (x - 1/8.d0 - 1/x + 5/xx/8)/3*b**2 -
     .        (x - 5/2.d0 + 2/x/3 + 5/xx/6)/4*b - zeta3/2/xx +
     .        zeta2*(1/2.d0 + 1/x - 5/xx/4) + 13/12.d0 - 5/x/24
      endif

      return
      end


C   Real part of two-point function B0 (mu is the t'Hooft scale), 
C   following up to an overall factor -16*pi2 the appendix of
C   G. Degrassi and A. Sirlin, PRD 46 (1992) 3104.

      double precision function b0(mu,x,y,z)
      
      implicit none
      double precision mu,m2,x,y,z,delta,cds,lambda,omega

         m2 = mu**2
      delta = (y - z)/x
        cds = (y + z)/x/2 - delta**2/4 - 1/4.d0

      if (y.eq.z) then
         b0 = 2*(1 - lambda(x/y)) - dlog(dabs(x/m2)) - dlog(dabs(y/x)) 
      else if (y.eq.0.d0) then
         b0 = 2 + (z/x - 1)*dlog(abs(1 - x/z)) - dlog(z/m2) 
      else if (z.eq.0.d0) then
         b0 = 2 + (y/x - 1)*dlog(abs(1 - x/y)) - dlog(y/m2) 
      else
         b0 = 2*(1 - omega(cds,delta)) - dlog(dabs(y/x))*(1 + delta)/2
     .              - dlog(dabs(x/m2)) - dlog(dabs(z/x))*(1 - delta)/2
      endif

      return
      end


      double precision function lambda(x)

C   A. Sirlin, NPB 332 (1990) 20 and 
C   G. Degrassi and A. Sirlin, PRD 46 (1992) 3104.

      implicit none
      double precision x,root

      if ((x.lt.4.d0).and.(x.gt.0)) then
         root = dsqrt(4/x - 1)
         lambda = root*datan(1/root)
      else
         root = dsqrt(dabs(4/x - 1))
         lambda = root/2*dlog(dabs((1 + root)/(1 - root)))
      endif

      return
      end


      double precision function omega(x,y)
      
C   G. Degrassi and A. Sirlin, PRD 46 (1992) 3104.

      implicit none
      double precision x,y,z,rx

      if (x.gt.0) then
            rx = dsqrt(x)
         omega = rx*(datan((1 + y)/2/rx) + datan((1 - y)/2/rx))
      else
            rx = dabs(dsqrt(x))
             z = x + y**2/4 - 1/4.d0
         omega = rx/2*dlog(dabs((z - rx)/(z + rx)))
      endif

      return
      end
