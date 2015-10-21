      subroutine anomagmntmu(amu)

C  computes a_mu = (g_mu - 2 - alpha/pi)/2.

      implicit none
      integer i,j
      dimension bbn(4)
      double precision cut,omega,weight,msrun,mcrun,mudeff,mseff,ln2
      double precision x,al,al3,arun,alfamu,alphas,as,as2,bbn,valdah
      double precision amu,a2tot,a2mu,a2e,a2tau,a2b,a2c,a2uds,a2had
      double precision a3tot,a3mvac,a3e,a3mu,a3tau,a3etau,a3elbl,a3mlbl
      double precision a3tlbl,a3muha,a3ehad,a3haha,a3hlbl,a4tot,a4mu,a4e
      double precision a4etau,a5estm,amuew,c(4,11),ag(2,11),c1,c2,cl2
      double precision rmemmu,rmuta2,rmuct2,rmuct4,rmumc2,rmumb2,rctmc2
      double precision cum0a2,cum0a3,cmatch,c3m,sin2tw,li41o2,s2,dgamma
      double precision lrmuct,lrmemu,lrmcct,lrmumc,lrmemc,lmzmmu,lmzmta
      double precision lnmzmb,lnmzmc,lmzmud,k1,kk1,k2,kk2,k2auds,k2buds
      double precision k2ac,k2bc,k2ab,k2bb

      include 'common.f'

      k1(x) = ( 25/ 4.d0 +  3*dlog(x))*x
      k2(x) = (291/10.d0 + 18*dlog(x))*x**2

      li41o2 = 0.5174790616738993863307582d0
          s2 = 4*cl2(pi1/3)/9/dsqrt(3.d0)
         ln2 = dlog(2.d0)

         cut =  1.800d0
      cmatch = 11/72.d0

      call coefhq(c,ag)

      mudeff = 176.d-3
       mseff = 305.d-3

       omega = 3.31849d0

C  1-loop result removed from definition!
C
C     a1 = 1/2.d0
C
C  old result and old implementation based on hep-ph/9805470:
C
C     a2had  =   63.430d-9 + 3.2d-6*(dahad3 - valdah)
C
C  The value 63.830 and the factor 1.006 below are from my average of 
C  the tau and e+e- based on results of hep-ph/0308213.

c     a2had  =   64.630d-9  ! PDG 2006 
      a2had  =   63.830d-9 

      a3muha = -  1.799d-9*1.006d0 
      a3ehad =    0.969d-9*1.006d0
      a3haha =    0.027d-9
      a3hlbl =    1.370d-9 ! +- 0.150d-9 (hep-ph/0605052)

      a3e    =    1.920455130d0 ! (33)   (hep-ph/0411168)
      a3tau  = -  0.00178233d0  ! (48)   (hep-ph/0411168)
      a3elbl =   20.94792489d0  ! (16)   (hep-ph/0411168)

      a4mu   = -   1.7260d0   ! +- 5.0d-3 Nucl.Phys.Proc.Suppl.144,206,2005  
      a4e    =   132.6823d0   ! +- 7.2d-3 (hep-ph/0402206)
      a4etau =     0.037594d0 ! +- 8.3d-5 (hep-ph/0402206)

C  5-loop estimate: compare with a5estm = 658.d0 hep-ph/0507174 (Kataev)

      a5estm =   652.0000d0   ! +- 20.0d0 (hep-ph/0512330)

      al = alpha/pi1
      al3 = al**3
      as = alphas(mc)/pi1
      if (f4lqcd.eqv..false.) then
         as = as*(1 + cmatch*as**2)
      else
         as = as*(1 + cmatch*as**2 + c3m(3)*as**3)
      endif
      as  = arun(3,mc,cut,as*pi1)/pi1
      as2 = as*as

      rmemmu =  me/mmu
      rmuta2 = (mmu/mtau)**2
      rmuct2 = (mmu/cut )**2
      rmuct4 =     rmuct2**2
      rmumc2 = (mmu/mc  )**2
      rmumb2 = (mmu/mb  )**2
      rctmc2 = (cut/mcrun(cut))**2

      lrmuct = dlog(rmuct2)
      lrmemu = dlog(rmemmu)
      lrmcct = dlog(rctmc2)
      lrmumc = dlog(rmumc2)
      lrmemc = dlog(me/mc)
      lmzmmu = dlog(mz/mmu)
      lmzmta = dlog(mz/mtau)
      lnmzmb = dlog(mz/mb)
      lnmzmc = dlog(mz/mc)
      lmzmud = dlog(mz/mudeff)

      sin2tw = 1 - 1/ratzw2

C  2-loop QED contribution from
C  G. Li, R. Mendel, and M.A. Samuel, PRD 47, 1723 (1993);
C  A. Czarnecki and M. Skrzypek, hep-ph/9812394.

      a2mu  = 197/144.d0 + 3*zeta3/4 + zeta2/2 - 3*zeta2*ln2
      a2e   = - 25/36.d0 - dlog(rmemmu)/3 + pi2*rmemmu/4
     .      + (3 + 4*dlog(rmemmu))*rmemmu**2 - 5*pi2/4*rmemmu**3
      a2tau = rmuta2/45 + rmuta2**2*(9/140.d0 + dlog(rmuta2))/140

C  light quark (uds) contribution:

      kk1 = k1(rmuct2) - 3*rmuct2/2
      kk2 = k2(rmuct2) - 6*rmuct4

      a2uds = 2*rmuct2*((1 + kk1/2 + kk2/3)*(1 + as) 
     .       + (( 17/ 3.d0 - 9*zeta3/2)*kk1 + 27*rmuct2/32 
     .       +  (281/72.d0 - 3*zeta3  )*kk2 +  3*rmuct4/2)*as2)/9

      if (flgech.eqv..true.) then
         j = 1.d4*as
         if (j.lt.0) then
            do 20 i = 1, 4
               bbn(i) = 0.d0
 20         continue
         else if (j.ge.2000) then
            do 30 i = 1, 4
               bbn(i) = bn(i,2000)
 30         continue
         else
            do 40 i = 1, 4
               bbn(i) = bn(i,j)+(as-asgrid(j))*(bn(i,j+1)-bn(i,j))*1.d4
 40         continue
         endif         
         cum0a2 =   299/24.d0  -   9*zeta3
         cum0a3 = 58057/288.d0 - 779*zeta3/4 + 75*zeta5/2 
         a2uds = a2uds + 2*rmuct2*((bbn(1) - as) + cum0a2*bbn(2) 
     .                                           + cum0a3*bbn(3))/9
      else
         cum0a2 =   245/24.d0  -   9*zeta3
         cum0a3 = 43675/288.d0 - 617*zeta3/4 + 75*zeta5/2 - 81*zeta2/8
         a2uds = a2uds + 2*rmuct2*as2*(cum0a2 + cum0a3*as)/9
      endif

      k2auds = al3/9*rmuct2*( - 4*pi2 + 23*lrmuct/3 + 377/9.d0 
     .       + (677*lrmuct/3 + 19*lrmuct**2 + 23647/24.d0 
     .       - 111*pi2)*rmuct2/24)*(1 + as)

      k2buds = al3/9*rmuct2*( - 8*lrmemu/3 - 2/3.d0 - (7*lrmuct/3 
     .       + lrmuct**2 + 137/24.d0 - pi2/3 + 4*lrmuct*lrmemu 
     .       + 19*lrmemu/3)*rmuct2)*(1 + as)

C  strange quark mass effects to uds contribution:

      a2uds = a2uds + 2*rmuct2*(msrun(cut)/cut)**2*as/9

C  charm quark mass effects (gluon splitting) to uds contribution:

      a2uds = a2uds + (mmu/mcrun(cut))**2*((3503/75.d0 - 2*pi2/3 
     .       - 88*lrmcct/5 + 2*lrmcct**2)/3 
     .       + rctmc2*(1723/420.d0 - lrmcct)/28)/405*as2

C  charm quark effects:

      as  = alphas(mc)/pi1
      as2 = as**2

      c1 = c(1,1) + c(2,1)*as + (c(3,1) + 3*c(4,1))*as2
      c2 = 108/1225.d0 - 0.1943d0*as + 3*dlog(rmumc2)*(c(1,2) 
     .   + c(2,2)*as + (c(3,2) + 3*c(4,2))*as2)

      a2c = rmumc2/36*(c1 + rmumc2/4*c2)

      k2ac = al3/6*(c1*(223/54.d0 - 2*zeta2 + 23*lrmumc/36) 
     .     - 3542/2025.d0)*rmumc2

      k2bc = - al3/108*c1*(1 + 4*lrmemu)*rmumc2

C  bottom quark effects:

      as  = alphas(mb)/pi1
      as2 = as**2

      c1 = c(1,1) + c(2,1)*as + (c(3,1) + 4*c(4,1))*as2
      c2 = 108/1225.d0 - 0.1943d0*as + 3*dlog(rmumb2)*(c(1,2) 
     .   + c(2,2)*as + (c(3,2) + 4*c(4,2))*as2)

      a2b = rmumb2/144*(c1 + rmumb2/4*c2)

      k2ab = al3/24*(c1*(223/54.d0 - 2*zeta2 + 23*dlog(rmumb2)/36)
     .     - 3542/2025.d0)*rmumb2

      k2bb = - al3/432*c1*(1 + 4*lrmemu)*rmumb2

      a2tot = al**2*(a2mu + a2e + a2tau + a2uds + a2c + a2b) + a2had

C  3-loop contribution:

      a3mlbl = (5 - 12*zeta3 + 15*zeta5/2 + 931*zeta2 - 123*zeta4/2
     .       - 15*zeta2*zeta3 - 6*pi2*ln2*(36 + ln2) + 6*ln2**4)/9 
     .       + 16*li41o2

      a3mu = 28259/5184.d0 + 139*zeta3/18 - (215*zeta5 + 239*zeta4)/24
     .     + 17101*zeta2/135 + 83*zeta2*zeta3/12 - 298*pi2/9*ln2
     .     + 25*(24*li41o2 + ln2**4 - pi2*ln2**2)/18

      a3tlbl = rmuta2*(3*zeta3/2 - 19/16.d0)

      a3mvac = a3mu - a3mlbl 

      a3etau = - rmuta2*(1 + 4*dlog(rmemmu))/135 + 2*(me/mtau)**2/15

      a3muha = a3muha + k2auds + k2ac + k2ab
      a3ehad = a3ehad + k2buds + k2bc + k2bb

      a3tot = al**3*(a3mu + a3e + a3tau + a3elbl + a3tlbl + a3etau)
     .      + a3muha + a3ehad + a3haha + a3hlbl

C  4-loop contributions:

      a4tot = al**4*(a4mu + a4e + a4etau)

C  electroweak 1-loop contribution:

      amuew = gf*mmu**2/24/dsqrt(2.d0)/pi2*(5 + (1 - 4*sin2tw)**2)

C  fermionic electroweak 2-loop contributions following
C  A. Czarnecki, B. Krause, and W.J. Marciano, Phys. Rev.D52, 2619 (1995);
C  G. Degrassi and G.F. Giudice, Phys. Rev. D58 (1998) 053007;
C  we included a formula which extrapolates between small and large Higgs mass
C  expansions and which is exact for mh = mt.  
C  Leading logarithms of bosonic electroweak 2-loop effects following
C  A. Czarnecki, B. Krause, and W.J. Marciano, Phys. Rev. Lett.76, 3267 (1996).

      call alfahat(mmu,dgamma,alfamu)
      weight = dexp( - omega*ratth2)

      amuew = amuew - gf*mmu**2/3/dsqrt(2.d0)/pi2*alfamu/pi1*
     .        (41/8.d0 - pi2/3 + 7/sin2tw/16 + 15*rattw2/sin2tw/128
     .      + 3/sin2tw/16*dlog(rattw2) + dlog(rattz2)
     .      + 3*dlog((mudeff/mtau)**3*(mc**2/mz)**2/mseff/mb)/4*
     .        (1 + 2*(1 - 4*sin2tw)/27) - 2*(1 - 4*sin2tw)/27*
     .        dlog(mudeff**7*(mc**2/mz)**4/mtau**9/mseff/mb)
     .      + (13 + 6*dlog(ratth2))/9*(1 - weight)
     .      + ratth2*(3 + pi2/3 + (dlog(ratth2) + 1)**2)*weight
     .      + lmzmmu*(123 - 124*sin2tw + 248*sin2tw**2)/12)

C  leading logarithms of 3-loop effects from
C  G. Degrassi and G.F. Giudice, Phys. Rev. D58 (1998) 053007.

      amuew = amuew + gf*mmu**2/24/dsqrt(2.d0)/pi2*alfamu**2/pi2*
     .        ((2827*lmzmmu**2/2 - 298*lmzmta**2 - 7826*lnmzmb**2/81)/9
     .      + 35200*lnmzmc**2/729               + 2108*lmzmud**2/81 
     .      + 24*(lnmzmb + 3*lmzmta - 4*lnmzmc - 2*lmzmud)*lmzmmu
     .      - 128*lnmzmb*lnmzmc/243 - 179*(lnmzmb**2/3 + lmzmta**2 
     .      +          4*lnmzmc**2/3 +   2*lmzmud**2 + 2*lmzmmu**2)/9
     .      +  4*(dlog(  mb/mtau  )**2/2 + 2*dlog(    mb/ mc)**2/3
     .      +     dlog(  mb/mudeff)**2   +   dlog(    mb/mmu)**2)
     .      - 16*(dlog(  mc/mudeff)**2   +   dlog(    mc/mmu)**2)
     .      + 12*(dlog(mtau/mudeff)**2   +   dlog(  mtau/mmu)**2)
     .      +  8*(dlog(mtau/mc    )**2   -   dlog(mudeff/mmu)**2))

      amu = (a2tot + a3tot + a4tot + al**5*a5estm + amuew)*1.d9

      return
      end


      subroutine coefhq(c,ag)

C  coefficients of heavy quark vector current correlator from 
C  K.G. Chetyrkin, J.H. K"uhn, and M. Steinhauser, Nucl. Phys.B505, 40 (1997);
C  converted to ms-bar mass definition. 

      double precision c(4,11),ag(2,11)

      c(1,1)  =      16/15.d0
      c(1,2)  =      16/35.d0
      c(1,3)  =     256/945.d0
      c(1,4)  =     128/693.d0
      c(1,5)  =    2048/15015.d0
      c(1,6)  =    2048/19305.d0
      c(1,7)  =   65536/765765.d0
      c(1,8)  =   16384/230945.d0
      c(1,9)  =  524288/8729721.d0
      c(1,10) =  524288/10140585.d0
      c(1,11) = 8388608/185910725.d0
      
      c(2,1)  =             3104.d0/1215.d0
      c(2,2)  =            15728.d0/14175.d0
      c(2,3)  =           773056.d0/1488375.d0
      c(2,4)  =           855136.d0/4209975.d0
      c(2,5)  =        522492928.d0/49165491375.d0
      c(2,6)  = -     6346221824.d0/54784404675.d0
      c(2,7)  = -  3976085794816.d0/19558032468975.d0
      c(2,8)  = - 10259619316736.d0/38566816162875.d0
      c(2,9)  = - 0.3122d0
      c(2,10) = - 0.3470d0
      c(2,11) = - 0.3735d0

C  O(alpha_s^2) part of Table 1; n_f independent part (exact results):
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
      c(3,9)  =   c(3,8)    ! estimate
      c(3,10) =   c(3,8)    ! estimate
      c(3,11) =   c(3,8)    ! estimate

C  O(alpha_s^2) part of Table 1; n_f dependent part:

      c(4,1)  =                2414.d0/3645.d0
      c(4,2)  =              290179.d0/637875.d0
      c(4,3)  =            22340824.d0/52093125.d0
      c(4,4)  =          1690158224.d0/3978426375.d0
      c(4,5)  =     216772869480064.d0/511075282843125.d0
      c(4,6)  =     626879108177216.d0/1480658105151225.d0
      c(4,7)  = 1115100660148360448.d0/2642974717694936625d0
      c(4,8)  = 7438303706248665728.d0/17719870787378384625.d0
      c(4,9)  = 0.418d0     ! estimate
      c(4,10) = 0.415d0     ! estimate
      c(4,11) = 0.413d0     ! estimate

C  gluon condensate contributions:

      ag(1,1)  = -       32/105.d0
      ag(1,2)  = -       32/63.d0
      ag(1,3)  = -      512/693.d0
      ag(1,4)  = -     1280/1287.d0
      ag(1,5)  = -     8192/6435.d0
      ag(1,6)  = -    57344/36465.d0
      ag(1,7)  = -   262144/138567.d0
      ag(1,8)  = -    65536/29393.d0
      ag(1,9)  = -  5242880/2028117.d0
      ag(1,10) = - 11534336/3900225.d0
      ag(1,11) = - 16777216/5014575.d0

      ag(2,1)  =             32099.d0/12960.d0
      ag(2,2)  =                59.d0/56.d0
      ag(2,3)  = -           20579.d0/42525.d0
      ag(2,4)  = -       100360567.d0/47628000.d0
      ag(2,5)  = -       459884251.d0/121080960.d0
      ag(2,6)  = -     63441631703.d0/11442150720.d0
      ag(2,7)  = -     31858548737.d0/4341887550.d0
      ag(2,8)  = -    197618805581.d0/21549939840.d0
      ag(2,9)  = -  13467196352539.d0/1220106888000.d0
      ag(2,10) = - 110308921063994.d0/8527479553593.d0
      ag(2,11) = - 14.861d0

      return
      end
