      subroutine smrule(sumrlc,momntc,errorc,sumrlb,momntb,errorb,ccmbc)

C   For the zeroth moment the gluon condensate contribution is
C   currently neglected.

      implicit none
      double precision alfamc,alfamb,dgamma,alphas,mcrun,mbrun,md0,mbpm
      double precision md02,mbpm2,asb,asc,as2m,lc,lb,lerr,sumrl,error
      double precision logm,agg,aggc,aggb,agger,aggerc,aggerb
      double precision fqcd1,fqcdb,fqcd2,root,logr,ccmbc,sumrlc(0:8)
      double precision sumrlb(0:8),momntc(8),momntb(8),errorc(0:8)
      double precision errorb(0:8),M_c(76),M_b(6),widthc(76),widthb(6)
      double precision widerc(76),widerb(6),c(4,8),b(2,8),coefc1(0:10)
      double precision coefb1(0:10),coeff2(8),coefc3(0:8),coefb3(0:8)
      double precision rtc1,rtc2,rtc3,rtc4,rtc5,rtc6,rtc7,rtc8,rtc9
      double precision rtb1,rtb2,rtb3,rtb4,rtb5,rtb6,rtb7,rtb8,rtb9
      integer k,n
      logical flgbes
      include 'common.f'

      data M_c/  3.096916d0, 3.68609d0, 
     .    3.750d0, 3.760d0, 3.764d0, 3.768d0, 3.770d0, 3.772d0, 3.776d0,
     .    3.780d0, 3.790d0, 3.810d0, 3.850d0, 3.890d0, 3.930d0, 3.940d0,
     .    3.950d0, 3.960d0, 3.970d0, 3.980d0, 3.990d0, 4.000d0, 4.010d0,
     .    4.020d0, 4.027d0, 4.030d0, 4.033d0, 4.040d0, 4.050d0, 4.060d0,
     .    4.070d0, 4.080d0, 4.090d0, 4.100d0, 4.110d0, 4.120d0, 4.130d0,
     .    4.140d0, 4.150d0, 4.160d0, 4.170d0, 4.180d0, 4.190d0, 4.200d0,
     .    4.210d0, 4.220d0, 4.230d0, 4.240d0, 4.245d0, 4.250d0, 4.255d0,
     .    4.260d0, 4.265d0, 4.270d0, 4.280d0, 4.300d0, 4.320d0, 4.340d0,
     .    4.350d0, 4.360d0, 4.380d0, 4.390d0, 4.400d0, 4.410d0, 4.420d0,
     .    4.430d0, 4.440d0, 4.450d0, 4.460d0, 4.480d0, 4.500d0, 4.520d0,
     .    4.540d0, 4.560d0, 4.600d0, 4.800d0/

      data widthc/  5.55d-6, 2.43d-6,
     .      3.5d0,  3.7d0,   4.2d0,   4.7d0,  2.6d0,  2.6d0,  4.1d0, 
     .      7.3d0,  5.7d0,   4.3d0,   9.3d0, 16.1d0, 23.6d0,  7.0d0, 
     .      7.3d0,  5.5d0,  10.5d0,   8.9d0,  8.2d0,  9.2d0, 12.9d0,
     .     18.7d0, 11.7d0,   7.0d0,  10.4d0, 18.4d0, 20.0d0, 24.2d0,
     .     19.1d0, 20.1d0,  18.3d0,  17.4d0, 16.9d0, 18.8d0, 17.6d0,
     .     16.0d0, 19.8d0,  18.9d0,  18.9d0, 19.5d0, 17.8d0, 16.4d0,
     .      9.7d0, 13.9d0,   9.8d0,   7.6d0,  3.7d0,  2.4d0,  3.2d0,
     .      3.7d0,  4.0d0,   7.7d0,  12.7d0, 17.6d0, 14.6d0, 15.6d0,
     .     12.6d0, 18.6d0,  19.0d0,  12.4d0, 16.8d0, 15.6d0, 14.5d0,
     .     17.9d0, 16.2d0,  15.2d0,  21.5d0, 26.2d0, 25.2d0, 20.4d0,
     .     20.0d0, 41.7d0, 129.8d0, 143.3d0/

      data widerc/  0.14d-6, 0.05d-6,
     .      1.8d0,  1.0d0,   1.5d0,   1.2d0,  0.4d0,  0.9d0,  1.2d0,
     .      0.9d0,  1.8d0,   4.0d0,   5.9d0,  6.5d0,  4.8d0,  2.1d0,
     .      1.9d0,  1.8d0,   1.4d0,   1.8d0,  2.1d0,  1.7d0,  2.2d0,
     .      1.8d0,  1.1d0,   0.8d0,   1.2d0,  1.7d0,  2.4d0,  2.1d0,  
     .      2.3d0,  2.3d0,   2.0d0,   2.0d0,  2.1d0,  3.0d0,  1.8d0,
     .      1.9d0,  2.2d0,   1.7d0,   1.9d0,  2.0d0,  1.4d0,  1.8d0,
     .      2.0d0,  2.2d0,   1.6d0,   1.2d0,  0.7d0,  0.7d0,  0.7d0,
     .      0.7d0,  0.8d0,   1.2d0,   2.4d0,  2.6d0,  3.1d0,  3.1d0,
     .      1.6d0,  2.8d0,   2.9d0,   1.9d0,  2.1d0,  2.1d0,  1.8d0,
     .      2.1d0,  2.0d0,   1.8d0,   2.9d0,  4.3d0,  3.3d0,  3.3d0,
     .      4.0d0,  4.9d0,  18.7d0,  19.8d0/

      flgbes = .false.

      md0   = 1.8645d0
      mbpm  = 5.2790d0
      md02  = 4*md0**2
      mbpm2 = 4*mbpm**2

      call alfahat(mc,dgamma,alfamc)
      call alfahat(mb,dgamma,alfamb)

      asc = alphas(mc)/pi1
      asb = alphas(mb)/pi1

      agg  = 0.0d0
      aggc = 2*zeta2*agg/mc**4
      aggb = 2*zeta2*agg/mb**4
      
      agger  = 0.07d0
      aggerc = 2*zeta2*agger/mc**4
      aggerb = 2*zeta2*agger/mb**4

C     lerr = 1.47d0 
      lerr = 1.24d0

C  The moment coefficients are normalized as they appear in the paper by
C  Chetyrkin, K"uhn, and Steinhauser, NPB 505 (1997) 40. For the "canonical"
C  normalization, M_n := m^2n integral R(s)/s^(n+1) ds, with asymptotically
C  R(s) = 1, they have to be multiplied by 3/4 (which is done at the end). 
C  In the canonical normalization all coefficients are of order unity and
C  higher orders can be estimated by considering group theory factors. 

      c(1,1)  =      16/15.d0
      c(1,2)  =      16/35.d0
      c(1,3)  =     256/945.d0
      c(1,4)  =     128/693.d0
      c(1,5)  =    2048/15015.d0
      c(1,6)  =    2048/19305.d0
      c(1,7)  =   65536/765765.d0
      c(1,8)  =   16384/230945.d0
C     c(1,9)  =  524288/8729721.d0
C     c(1,10) =  524288/10140585.d0
C     c(1,11) = 8388608/185910725.d0
      
      c(2,1)  =             3104.d0/1215.d0
      c(2,2)  =            15728.d0/14175.d0
      c(2,3)  =           773056.d0/1488375.d0
      c(2,4)  =           855136.d0/4209975.d0
      c(2,5)  =        522492928.d0/49165491375.d0
      c(2,6)  = -     6346221824.d0/54784404675.d0
      c(2,7)  = -  3976085794816.d0/19558032468975.d0
      c(2,8)  = - 10259619316736.d0/38566816162875.d0
C     c(2,9)  = - 0.3122d0
C     c(2,10) = - 0.3470d0
C     c(2,11) = - 0.3735d0

C   exact results:
C
C     c(3,1)  = - 1645049/34992 + 1537079/38880*zeta3
C     c(3,2)  = - 845815487803/2939328000 + 139672789/580608*zeta3
C     c(3,3)  = - 6185275180966013/3840721920000
C               + 311213994701/232243200*zeta3
C     c(3,4)  = - 369474963478061143/42247941120000 
C               + 6195053931629/851558400*zeta3
C     c(3,5)  = - 628032870023659409086469/13702021255802880000
C               + 106370602648560077/2789705318400*zeta3
C     c(3,6)  = - 674888245528059346325320577251/2898154099431828357120000
C               + 42569440523024104021/219742942003200*zeta3
C     c(3,7)  = - 16365348506250339255627072176208589/
C                 14189362470818231636459520000
C               + 1917072315222405358071049/1998028396088524800*zeta3
C     c(3,8)  = - 30514596739199641092015827094073861932059/
C                 5453965252908403694105945702400000 
C               + 2473736474011035386120449171/531475553359547596800*zeta3

      c(3,1)  =   0.5098816982876492d0 
      c(3,2)  =   1.4122711771465944d0 
      c(3,3)  =   0.3522281361474491d0
      c(3,4)  = - 0.4789375595174697d0
      c(3,5)  = - 0.9943079101012512d0
      c(3,6)  = - 1.2631281086183005d0
      c(3,7)  = - 1.3516966474732030d0
      c(3,8)  = - 1.3089622200501027d0

      c(4,1)  =                2414.d0/3645.d0
      c(4,2)  =              290179.d0/637875.d0
      c(4,3)  =            22340824.d0/52093125.d0
      c(4,4)  =          1690158224.d0/3978426375.d0
      c(4,5)  =     216772869480064.d0/511075282843125.d0
      c(4,6)  =     626879108177216.d0/1480658105151225.d0
      c(4,7)  = 1115100660148360448.d0/2642974717694936625d0
      c(4,8)  = 7438303706248665728.d0/17719870787378384625.d0

      b(1,1)  = -        32/105.d0
      b(1,2)  = -        32/63.d0
      b(1,3)  = -       512/693.d0
      b(1,4)  = -      1280/1287.d0
      b(1,5)  = -      8192/6435.d0
      b(1,6)  = -     57344/36465.d0
      b(1,7)  = -    262144/138567.d0
      b(1,8)  = -     65536/29393.d0
C     b(1,9)  = -   5242880/2028117.d0
C     b(1,10) = - 11534336/3900225.d0
C     b(1,11) = - 16777216/5014575.d0

      b(2,1)  =             32099.d0/12960.d0
      b(2,2)  =                59.d0/56.d0
      b(2,3)  = -           20579.d0/42525.d0
      b(2,4)  = -       100360567.d0/47628000.d0
      b(2,5)  = -       459884251.d0/121080960.d0
      b(2,6)  = -     63441631703.d0/11442150720.d0
      b(2,7)  = -     31858548737.d0/4341887550.d0
      b(2,8)  = -    197618805581.d0/21549939840.d0
C     b(2,9)  = -  13467196352539.d0/1220106888000.d0
C     b(2,10) = - 110308921063994.d0/8527479553593.d0
C     b(2,11) = - 14.861d0

      if (flgbes.eqv..false.) then
         M_c(1)  =  3.096916d0
         M_c(2)  =  3.68609d0
         M_c(3)  =  3.76990d0
         M_c(4)  =  4.04000d0
         M_c(5)  =  4.15900d0
         M_c(6)  =  4.41500d0

C        widthc(1) = 5.26d-6
         widthc(1) = 5.55d-6
C        widthc(2) = 2.19d-6
         widthc(2) = 2.43d-6
         widthc(3) = 0.26d-6
C        widthc(3) = 0.204d-6  ! CLEO 2005 (not yet included)
         widthc(4) = 0.75d-6
         widthc(5) = 0.77d-6
         widthc(6) = 0.47d-6
      
C        widerc(1) = 0.37d-6
         widerc(1) = 0.14d-6
C        widerc(2) = 0.15d-6
         widerc(2) = 0.05d-6
         widerc(3) = 0.04d-6
C        widerc(3) = 0.034d-6  ! CLEO 2005 (not yet included)
         widerc(4) = 0.15d-6
         widerc(5) = 0.23d-6
         widerc(6) = 0.10d-6
      endif

      M_b(1)  =  9.46030d0 
      M_b(2)  = 10.02326d0
      M_b(3)  = 10.35520d0
      M_b(4)  = 10.57930d0
      M_b(5)  = 10.86500d0
      M_b(6)  = 11.01900d0

      widthb(1) = 1.340d-6
      widthb(2) = 0.612d-6
      widthb(3) = 0.443d-6
      widthb(4) = 0.248d-6
      widthb(5) = 0.310d-6
      widthb(6) = 0.130d-6
      
      widerb(1) = 0.018d-6
      widerb(2) = 0.011d-6
      widerb(3) = 0.008d-6
      widerb(4) = 0.031d-6
      widerb(5) = 0.070d-6
      widerb(6) = 0.030d-6

      rtc1 = (mcrun(2*md0)/md0)**2
      rtc2 = rtc1*rtc1
      rtc3 = rtc2*rtc1
      rtc4 = rtc3*rtc1
      rtc5 = rtc4*rtc1
      rtc6 = rtc5*rtc1
      rtc7 = rtc6*rtc1
      rtc8 = rtc7*rtc1
      rtc9 = rtc8*rtc1
      
      root = dsqrt(1 - rtc1)
      logr = rtc1*dlog((1 + root)/(1 - root))/root/2
      
      coefc1(0) = dlog(4/rtc1) + dlog(2/rtc1*(1 - root) - 1)/root
      coefc1(1) = (1 - logr)/root**2
      coefc1(2) = (rtc1 + 2 + (rtc1 - 4)*logr)/4/root**4
      coefc1(3) = ( - 3*rtc2 + 10*rtc1 + 8 - 3*(rtc2 - 4*rtc1 + 8)*logr)
     .     /24/root**6
      coefc1(4) = (15*rtc3 - 62*rtc2 + 104*rtc1 + 48 + 3*(5*rtc3 
     .     - 24*rtc2 + 48*rtc1 - 64)*logr)/192/root**8
      coefc1(5) = ( - 105*rtc4 + 530*rtc3 - 1096*rtc2 + 1232*rtc1 + 384 
     .     - 15*(7*rtc4 - 40*rtc3 + 96*rtc2 - 128*rtc1 + 128)*logr)
     .     /1920/root**10
      coefc1(6) = (315*rtc5 - 1890*rtc4 + 4768*rtc3 - 6576*rtc2 
     .     + 5568*rtc1 + 1280 + 15*(21*rtc5 - 140*rtc4 + 400*rtc3
     .     - 640*rtc2 + 640*rtc1 - 512)*logr)/7680/root**12
      coefc1(7) = ( - 1155*rtc6 + 8050*rtc5 - 24136*rtc4 + 40576*rtc3
     .     - 41984*rtc2 + 28544*rtc1 + 5120 - 35*(33*rtc6 - 252*rtc5
     .     + 840*rtc4 - 1600*rtc3 + 1920*rtc2 - 1536*rtc1 + 1024)*logr)
     .     /35840/root**14
      coefc1(8) = (15015*rtc7 - 119350*rtc6 + 415688*rtc5 - 830448*rtc4 
     .     + 1046656*rtc3 - 866560*rtc2 + 492544*rtc1 + 71680 + 35*
     .     (429*rtc7 - 3696*rtc6 + 14112*rtc5 - 31360*rtc4 + 44800*rtc3 
     .     - 43008*rtc2 + 28672*rtc1 - 16384)*logr)/573440/root**16

      coefc3(0) = 2*( - 1 + logr/rtc1)/root**2

      rtb1 = (mbrun(2*mbpm)/mbpm)**2
      rtb2 = rtb1*rtb1
      rtb3 = rtb2*rtb1
      rtb4 = rtb3*rtb1
      rtb5 = rtb4*rtb1
      rtb6 = rtb5*rtb1
      rtb7 = rtb6*rtb1
      rtb8 = rtb7*rtb1
      rtb9 = rtb8*rtb1
      
      root = dsqrt(1 - rtb1)
      logr = rtb1*dlog((1 + root)/(1 - root))/root/2
      
      coefb1(0) = dlog(4/rtb1) + dlog(2/rtb1*(1 - root) - 1)/root
      coefb1(1) = (1 - logr)/root**2
      coefb1(2) = (rtb1 + 2 + (rtb1 - 4)*logr)/4/root**4
      coefb1(3) = ( - 3*rtb2 + 10*rtb1 + 8 - 3*(rtb2 - 4*rtb1 + 8)*logr)
     .     /24/root**6
      coefb1(4) = (15*rtb3 - 62*rtb2 + 104*rtb1 + 48 + 3*(5*rtb3 
     .     - 24*rtb2 + 48*rtb1 - 64)*logr)/192/root**8
      coefb1(5) = ( - 105*rtb4 + 530*rtb3 - 1096*rtb2 + 1232*rtb1 + 384 
     .     - 15*(7*rtb4 - 40*rtb3 + 96*rtb2 - 128*rtb1 + 128)*logr)
     .     /1920/root**10
      coefb1(6) = (315*rtb5 - 1890*rtb4 + 4768*rtb3 - 6576*rtb2
     .     + 5568*rtb1 + 1280 + 15*(21*rtb5 - 140*rtb4 + 400*rtb3
     .     - 640*rtb2 + 640*rtb1 - 512)*logr)/7680/root**12
      coefb1(7) = ( - 1155*rtb6 + 8050*rtb5 - 24136*rtb4 + 40576*rtb3
     .     - 41984*rtb2 + 28544*rtb1 + 5120 - 35*(33*rtb6 - 252*rtb5
     .     + 840*rtb4 - 1600*rtb3 + 1920*rtb2 - 1536*rtb1 + 1024)*logr)
     .     /35840/root**14
      coefb1(8) = (15015*rtb7 - 119350*rtb6 + 415688*rtb5 - 830448*rtb4 
     .     + 1046656*rtb3 - 866560*rtb2 + 492544*rtb1 + 71680 + 35*
     .     (429*rtb7 - 3696*rtb6 + 14112*rtb5 - 31360*rtb4 + 44800*rtb3 
     .     - 43008*rtb2 + 28672*rtb1 - 16384)*logr)/573440/root**16

      coefb3(0) = 2*( - 1 + logr/rtb1)/root**2

      do 5 n = 1, 8
         coeff2(n) = 1.d0/n**2
         coefc3(n) = (coefc3(n-1) - coefc1(n))/(1 - rtc1)
         coefb3(n) = (coefb3(n-1) - coefb1(n))/(1 - rtb1)
 5    continue
      
      do 100 n = 0, 8
         
C  charm quark:
         
         as2m  = alphas(2*md0)/pi1
         logm  = 2*dlog(2*md0/mc)
         fqcd1 = 1 + as2m + as2m**2*(277/24.d0 - 25*zeta3/3)
         fqcd2 =       - 25*as2m**2/12

         sumrl = 0.d0
         error = 0.d0
         do 10 k = 1, 2
            if ((flgbes.eqv..true.).and.(k.ge.3)) then
               sumrl = sumrl +  
     .              widthc(k)/M_c(k)**(2*n + 1)*alpha**2/pi1/4500.d0
               error = error + 
     .             (widerc(k)/M_c(k)**(2*n + 1)*alpha**2/pi1/4500.d0)**2
            else
               sumrl = sumrl +  widthc(k)/M_c(k)**(2*n + 1)
               error = error + (widerc(k)/M_c(k)**(2*n + 1))**2
            endif
 10      continue
         
         sumrl =  9*pi1/alfamc**2*sumrl
         error = (9*pi1/alfamc**2)**2*error

         if (n.eq.0) then
            lc = (((fqcd1 - as2m)*logm - 12*dlog(1 - 25*as2m/12*logm)/25 
     .           - 5/3.d0 + asc*(4*zeta3 - 7/2.d0) + asc**2
     .           *((18213*zeta3 - 17471)/432 - 25*zeta5/3 + 25*zeta2/12) 
     .           - 3*sumrl/4)/fqcd1 - coefc1(0))*2/rtc1/coefc3(0)

            sumrlc(0) = sumrl 
     .           + 4*(fqcd1*(coefc1(0) + rtc1*coefc3(0)/2*lc) - ((fqcd1 
     .           - as2m)*logm - 12*dlog(1 - 25*as2m/12*logm)/25))/3

            errorc(0) = dsqrt(error + 256.d0*asc**6
     .           + (2*lerr/3*fqcd1*rtc1*coefc3(0))**2)
         else
            sumrlc(n) = sumrl + 4*(fqcd2*coeff2(n) + fqcd1*(coefc1(n) 
     .           + rtc1*coefc3(n)/2*lc))/md02**n/3

            momntc(n) = (c(1,n) + c(2,n)*asc + asc**2*(c(3,n) +3*c(4,n)) 
     .           + aggc*b(1,n)*(1 + b(2,n)*asc))/(2*mc)**(2*n)

            errorc(n) = dsqrt(error + (256.d0*asc**6
     .           + (aggerc*b(1,n)*(1 + b(2,n)*asc))**2)/(2*mc)**(4*n)
     .           + (2*lerr/3*fqcd1*rtc1*coefc3(n)/md02**n)**2)

            sumrlc(n) = sumrlc(n)*10.d0**n
            momntc(n) = momntc(n)*10.d0**n
            errorc(n) = errorc(n)*10.d0**n
         endif

C  bottom quark:

         as2m  = alphas(2*mbpm)/pi1
         logm  = 2*dlog(2*mbpm/mb)
         fqcdb = 1 + as2m + as2m**2*(85/8.d0 - 23*zeta3/3)
         fqcd2 =       - 23*as2m**2/12

         sumrl = 0.d0
         error = 0.d0
         do 11 k = 1, 3
            sumrl = sumrl +  widthb(k)/M_b(k)**(2*n + 1)
            error = error + (widerb(k)/M_b(k)**(2*n + 1))**2
 11      continue

         sumrl =  9*pi1/alfamb**2*sumrl
         error = (9*pi1/alfamb**2)**2*error

         if (n.eq.0) then
            lb = (((fqcdb - as2m)*logm - 12*dlog(1 - 23*as2m/12*logm)/23 
     .           - 5/3.d0 + asb*(4*zeta3 - 7/2.d0) + asb**2
     .           *((17301*zeta3 - 16117)/432 - 25*zeta5/3 + 23*zeta2/12) 
     .           - 3*sumrl)/fqcdb - coefb1(0))*2/rtb1/coefb3(0)

            sumrlb(0) = (sumrl 
     .           + (fqcdb*(coefb1(0) + rtb1*coefb3(0)/2*lb) - ((fqcdb 
     .           - as2m)*logm - 12*dlog(1 -23*as2m/12*logm)/23))/3)*1.d1

            errorb(0) = dsqrt(error + 16.d0*asb**6
     .           + (lerr/6*fqcdb*rtb1*coefb3(0))**2)
         else
            sumrlb(n) = sumrl + (fqcd2*coeff2(n) + fqcdb*(coefb1(n)
     .           + rtb1*coefb3(n)/2*lb))/mbpm2**n/3
 
            momntb(n) = (c(1,n) + c(2,n)*asb + (c(3,n) +4*c(4,n))*asb**2 
     .           + aggb*b(1,n)*(1 + b(2,n)*asb))/4/(2*mb)**(2*n)

            errorb(n) = dsqrt(error + (16.d0*asb**6
     .           + (aggerb*b(1,n)*(1 + b(2,n)*asb)/4)**2)/(2*mb)**(4*n)
     .           + (lerr/6*fqcdb*rtb1*coefb3(n)/mbpm2**n)**2)

            sumrlb(n) = sumrlb(n)*10.d0**(2*n+1)
            momntb(n) = momntb(n)*10.d0**(2*n+1)
            errorb(n) = errorb(n)*10.d0**(2*n+1)
         endif

 100  continue

      ccmbc = (c(3,2) + 3*c(4,2))/(c(3,6) + 4*c(4,6))
      if (dabs(ccmbc).gt.1) ccmbc = 1/ccmbc

      ccmbc = ((aggerc*b(1,2)*(1 + b(2,2)*asc)*aggerb*b(1,6)*(1 
     .     + b(2,6)*asb) + ccmbc*256*(asc*asb)**3)/(2*mc)**4/(2*mb)**12
     .     + (2*lerr/3)**2*fqcd1*fqcdb*rtc1*rtb1*coefc3(2)*coefb3(6)
     .     /md02**2/mbpm2**6)/errorc(2)/errorb(6)/4*1.d15

      return
      end
