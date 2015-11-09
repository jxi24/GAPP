      subroutine z0pole(gamma,sigmah,R,sin2te,alr,afb,fmin,fmax)

      implicit none
      integer fmin,fmax,f
      double precision gamma(0:11),gammav(0:9),gammaa(0:9),sigmah,R(0:9)
      double precision sin2te(0:9),alr(0:10),afb(0:10)
      double precision cum0a2,cum0a3,cum0a4,cum0a5,cum0a6
      double precision fbeta1,fbeta2,fbeta3,c1,c2,c3,cpi2,cpi4
      double precision cum2a3,cumta2,cumta3,cal2
      double precision        cvm2a1,cvm2a2,cvm2a3,cvm4a0,cvm4a1,cvm4a2
      double precision cam2a0,cam2a1,cam2a2,cam2a3,cam4a0,cam4a1,cam4a2
      double precision csaa2,csaa2t,csaa3t,csm2a2,csva3,csva3t
      double precision rtz2,rtz4,rtz6,sumv,sumaef,sinthe
      double precision mtauz2,mcz2,mbz2,mbz4,msrun,mcrun,mbrun,mtrun
      double precision al,az,az2,az3,alphas,lntz2,lntz22,qf234,hbarc2
      complex*16 rho(0:9),kappa(0:9),veff(0:9),aeff(0:9)

C Kai028 ---------------------------------------
      double precision gammaFirstSecond, gammaThird
C ----------------------------------------------

      include 'common.f'

      hbarc2 = 0.38937966d6

       al = alphat/pi1
       az = alphas(mz)/pi1
      az2 = az**2
      az3 = az**3

      mtauz2 =      (mtau/mz)**2
        mcz2 = (mcrun(mz)/mz)**2
        mbz2 = (mbrun(mz)/mz)**2
        mbz4 =           mbz2**2
        mcz2 = mcz2*(1 + 2/3.d0*al*dlog(mcz2))
        mbz2 = mbz2*(1 + 1/6.d0*al*dlog(mbz2))

        rtz2 = (mtrun(mz)/mz)**2
        rtz4 =           rtz2**2
        rtz6 =           rtz2**3

      lntz2  = dlog(rtz2)
      lntz22 = dlog(rtz2)**2

      call rhof(rho,0,9)
      call kappaf(kappa,0,9)

      do 10 f = 0, 9
         veff(f) = cdsqrt(rho(f)/sinhat/coshat)/2
     .           * (i3(f) - 2*kappa(f)*q(f)*sinhat)
         aeff(f) = cdsqrt(rho(f)/sinhat/coshat)/2*i3(f)
         if (flagzp.ne.0) then
             sinthe = ratg21*sinth/2/dsqrt(sinhat*coshat)
            veff(f) = costh*veff(f) + sinthe*v2(f)
            aeff(f) = costh*aeff(f) + sinthe*a2(f)
         endif
         if (f.ne.0) sin2te(f) = (1-dreal(veff(f)/aeff(f)))/dabs(q(f))/4
 10   continue

      cal2   =       110/9.d0 -     80*zeta3/9
      cpi2   =                -    529*zeta2/72
      cpi4   =                  279841*zeta4/1152
      cum0a2 =        85/8.d0 -     23*zeta3/3
      cum0a3 = 372841/2592.d0 -  15961*zeta3/108  + 575*zeta5/18
      cum0a4 = (419585161 - 485655120*zeta3 + 79587360*zeta3**2  
     .     + 68382360*zeta5 - 19822320*zeta7)/186624
      cum0a5 = 0.d0
      cum0a6 = 0.d0

      if (flgech.eqv..true.) then
             c1 = 12*fbeta1(5)/23
             c2 = 12*fbeta2(5)/23
             c3 = 12*fbeta3(5)/23
         cum0a5 = cpi2*(6*cum0a3 + 7*c1*cum0a2 + 3*c1**2/2 + 3*c2) 
     .          + cpi4
         cum0a6 = cpi2*(27*c1/2*cum0a3 + 4*c1**2*cum0a2 + 7*c1*c2/2 
     .          +          8*c2*cum0a2 + 7*c3/2)
     .          + cpi4*(5*cum0a2 + 77*c1/12)
C  old effective charge or PMS optimization from Kataev and Starshenko:
C
C        cum0a4 = cum0a2*(9769/6624.d0 - 29*cum0a2/46 - 2*cum0a2**2 
C    .          + 3*cum0a3) - 3335*zeta2/144
      endif

      cum0a3 = cum0a3 + cpi2
      cum0a4 = cum0a4 + cpi2*(3*cum0a2 + 5*c1/2)

      cum2a3 = -     560/9.d0 +  140*zeta3/3
C     cum4a2 =        13/3.d0 -    4*zeta3                          (neglected)
      cumta2 =        (44/675.d0 + 2*lntz2/135  )/rtz2
     .       -  (1303/1058400.d0 +   lntz2/2520 )/rtz4
     .       + (1643/26790750.d0 +   lntz2/42525)/rtz6

      cumta3 = 
     . - (1.73743d-1 +  185827*lntz2/8.748d5   + 181*lntz22/4860  )/rtz2
     . - (7.52181d-3 -   19223*lntz2/3.26592d7 - 139*lntz22/362880)/rtz4
     . + (5.04106d-4 +18671173*lntz2/1.5431472d11
     .                                       + 481*lntz22/1.5309d7)/rtz6
C
C  exact coefficients extracted from 
C  S. A. Larin, T. van Ritbergen and J. A. M. Vermaseren, NIKHEF-H/94-30:
C
C  -        486587/5248800         +  49*zeta2/2430    -   1847*zeta3/19440
C  -     391913989/8001504000      - 163*zeta2/181440  + 103693*zeta3/2903040
C    4383502001359/194436547200000 + 299*zeta2/7654500 -   4943*zeta3/268800
                                                                      
      cvm2a1 =          12.d0
      cvm2a2 =       629/6.d0 
      cvm2a3 =    89893/54.d0 -  1645*zeta2/6 +   820*zeta3/27 
     .                                        - 36575*zeta5/54

      cvm4a0 = -         6.d0
C     cvm4a1 = -        22.d0                                       (neglected)
C     cvm4a2 = -   7657/36.d0 +   142*zeta2   +   296*zeta3/3       (neglected)
C     cvlna2 = -      23/6.d0                                       (neglected)

      cam2a0 = -         6.d0
      cam2a1 = -        22.d0
      cam2a2 = -    2237/8.d0 +    47*zeta2   +    97*zeta3
      cam2a3 = - 25024465/7776.d0             + 15515*zeta2/18 
     .       + 27545*zeta3/12 +    25*zeta4   -   995*zeta5

      cam4a0 =           6.d0
C     cam4a1 =          10.d0                                       (neglected)
C     cam4a2 =       993/4.d0 -   142*zeta2   -   580*zeta3/3       (neglected)
C     calna2 =       161/6.d0                                       (neglected)

      csaa2  = -      17/2.d0
      csaa2t =         3/4.d0 +    3*lntz2   
     .       -      7/rtz2/27 -    7/rtz4/400 - 104/rtz6/55125
 
      csaa3t = -   5651/72.d0 +    23*zeta2/2 +   3*zeta3
     .               -     31*lntz2/6        + 23*lntz22/4
     . - (3.77554d-1 +   6961*lntz2/8.1d3    +    lntz22/30  )/rtz2
     . - (9.39672d-3 + 115421*lntz2/3.1752d6 +    lntz22/280 )/rtz4 
     . - (2.66481d-3 + 135059*lntz2/9.5256d7 +    lntz22/1890)/rtz6
C
C  exact coefficients extracted from 
C  S. A. Larin, T. van Ritbergen and J. A. M. Vermaseren, NIKHEF-H/94-30:
C
C  -        7351031/1944000       + zeta2/15  +       3157*zeta3/1152
C  -    45837722011/8534937600    + zeta2/140 +    4921619*zeta3/1105920 
C  - 35954276563439/4389396480000 + zeta2/945 + 2343291463*zeta3/344064000
                                                                      
      csm2a2 =       18.d0 +  6*lntz2 
     .         - (80/81.d0 +  5*lntz2/27)/rtz2

      csva3  =                 55/72.d0 -    5*zeta3/3 
      csva3t = -            (185/288.d0 -   65*zeta3/108  )/rtz2
     .         -      (110401/233280.d0 -   29*zeta3/72   )/rtz4
     .         - (45317647/114307200.d0 - 5009*zeta3/15120)/rtz6

C   Universal corrections (QED, QCD, mixed and massive):

      do 20 f = 0, 9
         if (f.eq.6) then
            gammav(6) = 0.d0
            gammaa(6) = 0.d0
             gamma(6) = 0.d0
               alr(6) = 0.d0
               afb(6) = 0.d0
         else
            qf234 = q(f)**2*3/4
            gammav(f) = nc(f)/3*alphat*cdabs(veff(f)**2)*mz
            gammaa(f) = nc(f)/3*alphat*cdabs(aeff(f)**2)*mz
C            gamma(f) = (gammav(f) + gammaa(f))*(1 + al*qf234)

C Kai012 ---------------------------------------
C     Corrections to the Z-decay width 
C ----------------------------------------------

            gamma(f) = ( gammav(f) + gammaa(f) )
           
      if( modtype.eq.1 ) then

            gamma(f) = gamma(f) +
     -      (alphat*mz*nc(f)*(fits2b*
     -       ((-kkcc + kkss)*i3(f)**2 - 2*kkss**2*i3(f)*q(f) + 
     -         2*kkss**2*(kkcc + kkss)*q(f)**2) + 
     -      fitcph**2*((-kkcc + kkss)*i3(f)**2 + 
     -         2*(kkcc + (-1 + kkcc)*kkss)*i3(f)*q(f) - 
     -         2*kkss*((-2 + kkss)*kkss + kkcc*(2 + kkss))*
     -          q(f)**2) + 2*fitcph*(kkcc - kkss)*
     -       (i3(f)*(i3(f) - q(f)) + 
     -         kkss*q(f)*(-i3(f) - i3r(f) + 2*q(f)))))/
     -  (6.*fitx*kkcc*(kkcc - kkss)*kkss)

      endif

      if( modtype.eq.2 ) then

            gamma(f) = gamma(f) +
     -     (alphat*mz*nc(f)*(2*fits2b*
     -       ((-kkcc + kkss)*i3(f)**2 - 2*kkss**2*i3(f)*q(f) + 
     -         2*kkss**2*(kkcc + kkss)*q(f)**2) + 
     -      fitcph**2*((-kkcc + kkss)*i3(f)**2 + 
     -         2*(kkcc + (-1 + kkcc)*kkss)*i3(f)*q(f) - 
     -         2*kkss*((-2 + kkss)*kkss + kkcc*(2 + kkss))*
     -          q(f)**2) + 
     -      2*fitcph*(kkcc - kkss)*
     -       (i3(f)*(i3(f) - q(f)) + 
     -         kkss*q(f)*(-i3(f) - i3r(f) + 2*q(f)))))/
     -  (24.*fitx*kkcc*(kkcc - kkss)*kkss)

      endif

      if( modtype.eq.3 ) then

         if(f.le.3) then

            gamma(f) = gamma(f) + 
     -     (alphat*fitsph**2*mz*nc(f)*
     -    ((-kkcc + kkss)*i3(f)**2 + 2*kkcc*kkss*i3(f)*q(f) - 
     -      2*kkss**2*(kkcc + kkss)*q(f)**2))/
     -  (6.*fitx*kkcc*(kkcc - kkss)*kkss)

         else

            gamma(f) = gamma(f) +
     -     (alphat*fitsph*mz*nc(f)*
     -    ((-2 + fitsph)*(-kkcc + kkss)*i3(f)**2 + 
     -      2*kkss*((-1 + fitsph)*kkcc + kkss)*i3(f)*q(f) - 
     -      2*fitsph*kkss**2*(kkcc + kkss)*q(f)**2))/
     -  (6.*fitx*kkcc*(kkcc - kkss)*kkss)


         endif

      if( modtype.eq.4 ) then

         gammaFirstSecond = (alphat*fitsph*mz**((-2 + fitsph)*
     -     (-kkcc + kkss)*i3(f)**2 + 
     -     2*kkss*((-1 + fitsph)*kkcc + kkss)*i3(f)*q(f) - 
     -     2*fitsph*kkss**2*(kkcc + kkss)*q(f)**2))/
     -     (6.*fitx*kkcc*(kkcc - kkss)*kkss)

         gammaThird = (alphat*fitsph**2*mz*
     -     ((-kkcc + kkss)*i3(f)**2 + 2*kkcc*kkss*i3(f)*q(f) - 
     -     2*kkss**2*(kkcc + kkss)*q(f)**2))/
     -     (6.*fitx*kkcc*(kkcc - kkss)*kkss)

         if(f.eq.0) then

            gamma(f) = gamma(f) +
     -      2.0d0*gammaFirstSecond + gammaThird

         else

         if(f.eq.3.or.f.eq.9) then

            gamma(f) = gamma(f) + nc(f)*gammaThird

         else

            gamma(f) = gamma(f) + nc(f)*gammaFirstSecond

         endif

         endif

      endif

      endif

      gamma(f)=gamma(f)*(1 + al*qf234)

C Kai012 ---------------------------------------

             if (fa2mt0.eqv..true.) gamma(f) =  gamma(f) 
     .            - (gammav(f) + gammaa(f))*al**2*qf234*(qf234/6 + cal2)
             alr(f) = (1 - 8*i3(f)*q(f)*sin2te(f))/
     .                (1 - 8*i3(f)*q(f)*sin2te(f)+8*(q(f)*sin2te(f))**2)

C           write(50,*) 'alr ', f
C           write(50,*) alr(f)
C Kai013 ---------------------------------------
C     Corrections to the LR asymmetries
C ----------------------------------------------

      if( modtype.eq.1 ) then

             alr(f) = alr(f) +     
     -   (4*kkss*q(f)*(-i3(f) + kkss*q(f))*
     -    ((-1 + kkss)*(fitcph + fitcph**2*(-1 + kkss) - 
     -         2*fitcph*kkss + fits2b*kkss)*i3(f)*q(f) + 
     -      fitcph*(-1 + 2*kkss)*i3r(f)*(-i3(f) + kkss*q(f))))/
     -  (fitx*(-1 + 2*kkss)*
     -    (i3(f)**2 - 2*kkss*i3(f)*q(f) + 2*kkss**2*q(f)**2)**2)

      endif

      if( modtype.eq.2 ) then

             alr(f) = alr(f) + 
     -  (kkss*q(f)*(-i3(f) + kkss*q(f))*
     -    ((-1 + kkss)*(fitcph + fitcph**2*(-1 + kkss) - 
     -         2*fitcph*kkss + 2*fits2b*kkss)*i3(f)*q(f) + 
     -      fitcph*(-1 + 2*kkss)*i3r(f)*(-i3(f) + kkss*q(f))))/
     -  (fitx*(-1 + 2*kkss)*
     -    (i3(f)**2 - 2*kkss*i3(f)*q(f) + 2*kkss**2*q(f)**2)**2)
     
      endif

      if( modtype.eq.3 ) then

         if(f.le.3) then

            alr(f) = alr(f) + 
     -     (4*fitsph**2*kkss**3*i3(f)*q(f)**2*
     -     (i3(f) - kkss*q(f)))/(fitx*(kkcc - kkss)*
     -     (i3(f)**2 - 2*kkss*i3(f)*q(f) + 2*kkss**2*q(f)**2)**2)

C           write(50,*) alr(f) 
C           write(50,*) '========='

         else

            alr(f) = alr(f) + 
     -     (4*fitsph*kkss**2*(1 + (-2 + fitsph)*kkss)*i3(f)*
     -     q(f)**2*(-i3(f) + kkss*q(f)))/(fitx*(-1 + 2*kkss)*
     -     (i3(f)**2 - 2*kkss*i3(f)*q(f) + 2*kkss**2*q(f)**2)**2)

C           write(50,*) alr(f) 
C           write(50,*) '========='

         endif

      if( modtype.eq.4 ) then

         if(f.le.1) then

         if(f.eq.3.or.f.eq.9) then

            alr(f) = alr(f) + 
     -   (4*fitsph**2*kkss**3*i3(f)*q(f)**2*
     -   (i3(f) - kkss*q(f)))/ (fitx*(kkcc - kkss)*
     -   (i3(f)**2 - 2*kkss*i3(f)*q(f) + 2*kkss**2*q(f)**2)**2)

         else

            alr(f) = alr(f) + 
     -   (4*fitsph*kkss**2*(1 + (-2 + fitsph)*kkss)*i3(f)*
     -   q(f)**2*(-i3(f) + kkss*q(f)))/ (fitx*(-1 + 2*kkss)*
     -   (i3(f)**2 - 2*kkss*i3(f)*q(f) + 2*kkss**2*q(f)**2)**2)

         endif

         endif

      endif

      endif

C Kai013 ---------------------------------------

             if (f.eq.4) then
                afb(4) = 3/4.d0*alr(1)*(5 - 12*sin2te(4))/
     .               (5 - 12*sin2te(4) + 136*sin2te(4)**2/9)
             else
                afb(f) = 3/4.d0*alr(1)*alr(f)
             endif
             if (f.ge.4) gamma(f) =  gamma(f) + (gammav(f) + gammaa(f))*
     .                   (az + (cum0a2 + cumta2)*az2 - al*az/3*qf234
     .                 + (cum0a3 + cum2a3*(mcz2 + mbz2) + cumta3)*az3 
     .              + cum0a4*az2**2 + cum0a5*az2*az3 + cum0a6*az3*az3)
         endif
 20   continue

C   Massive non-singlet corrections:
      
      gamma(2) = gamma(2) + gammaa(2)*cam2a0*(mmu/mz)**2

      gamma(3) = gamma(3) + gammaa(3)*cam2a0*mtauz2
     .         * (1 - 2*al*(1 - 3/4.d0*dlog(mtauz2)))

      gamma(5) = gamma(5) + gammav(5)*mcz2*          cvm2a1*az
     .                    + gammaa(5)*mcz2*(cam2a0 + cam2a1*az)

      gamma(8) = gamma(8) + gammaa(8)*cam2a0*(msrun(mz)/mz)**2

      gamma(9) = gamma(9) 
     .         + gammav(9)*(mbz4*cvm4a0 + mbz2*
     .           (         cvm2a1*az + cvm2a2*az2 + cvm2a3*az3))
     .         + gammaa(9)*(mbz4*cam4a0 + mbz2*
     .           (cam2a0 + cam2a1*az + cam2a2*az2 + cam2a3*az3))

C   Singlet corrections:

      sumv   = - (1 + 4*sinhat/3)/4/dsqrt(sinhat*coshat)
      sumaef = 0.d0
      do 30 f = 4, 9
         if (f.ne.6) sumaef = sumaef + dreal(aeff(f))
 30   continue

      do 40 f = 4, 9
         if (f.ne.6) then
            gamma(f) = gamma(f) + az3*alphat/3*mz*v(f)*
     .                 (sumv*csva3 + v(6)*csva3t)
            gamma(f) = gamma(f) + az2*alphat/3*mz*dreal(aeff(f))*
     .                 (sumaef*csaa2 + dreal(aeff(6))*csaa2t) 
            gamma(f) = gamma(f) - az3*alphat/3*mz*a(f)*a(6)*csaa3t
         endif
 40   continue

      gamma(9) = gamma(9) - az2*alphat*mz*mbz2*a(9)*a(6)*csm2a2

      if (fobliq.eqv..true.) gamma(0) = gamma(0)*Zpar/3.d0
      
      gamma(10) = 0.d0
        alr(10) = 0.d0
      do 60 f = 4, 9
         gamma(f)  = 1.d3*gamma(f)
         gamma(10) = gamma(10) + gamma(f)
           alr(10) = alr(10) - q(f)/dabs(q(f))*gamma(f)*alr(f)
 60   continue

      gamma(11) = gamma(10)
      do 70 f = 0, 3
         gamma(f) = 1.d3*gamma(f)
         gamma(11)   = gamma(11) + gamma(f)
 70   continue

      sigmah = 12*pi1*gamma(1)*gamma(10)/(mz*gamma(11))**2*hbarc2

      do 80 f = 0, 3
         R(f) = gamma(10)/gamma(f)
 80   continue

      do 90 f = 4, 9
         R(f) = gamma(f)/gamma(10)
 90   continue

      alr(10) = alr(10)/gamma(10)
      afb(10) = 3/4.d0*alr(1)*alr(10)

      return
      end
