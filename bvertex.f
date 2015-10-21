      subroutine bvertex(dkbtop,taub)

C     special Zbb-vertex correction following 
C     J. Bernabeu, A. Pich and A Santamaria, NPB 363 (1991) 326.

      implicit none
      double precision x,y,z,lambda,alphas,at,gds,xasmt2
      double precision taub,taub2,fftj0,fftj1,gftj,xt,aa
      complex*16 r1,r2,r3,r4,f0a,f1a,f2a,f0b,f1b,f2b,spence,fds
      complex*16 i1a,i1a0,i1b,i1b0,i1cd,i1ef,i2a,i2b,i2cd,fb,dkbtop
      include 'ff.h'
      include 'common.f'

      x  = rattw2
      y  = ratzw2
      z  = 1 - x

      r1 = (1 + (0,1)*dsqrt(4/y*x - 1))/2
      r2 = (1 - (0,1)*dsqrt(4/y*x - 1))/2
      r3 = (1 + (0,1)*dsqrt(4/y   - 1))/2
      r4 = (1 - (0,1)*dsqrt(4/y   - 1))/2

      f0a = 2*(lambda(y/x) - 1) + dlog(x)
      f0b = 2*(lambda(y)   - 1)
      f1a = 2/y*(  spence(1/r1)          + spence(1/r2) 
     .           - spence(1/(r1 + x*r2)) - spence(1/(r2 + x*r1)) 
     .           + spence(x/(r1 + x*r2)) + spence(x/(r2 + x*r1)) 
     .           - dlog(x) * dlog(1 + y/z**2))
      f1b = 2/y*(  spence(1/r3)          + spence(1/r4)
     .           + spence(1/(r3 + x*r4)) + spence(1/(r4 + x*r3))
     .           - spence(x/(r3 + x*r4)) - spence(x/(r4 + x*r3)) 
     .           + dlog(x) * dlog(1 + x*y/z**2))
      f2a = -4/y*(f0a + 1 + x*dlog(x)/z + z/2*f1a)
      f2b = -4/y*(f0b + 1 + x*dlog(x)/z - z/2*f1b)

      i1a  = (3 - 4*sinhat)/6*(1.5d0 + f0a 
     .       + (y + z/2)*f2a - y*f1a) + 2/3.d0*x*sinhat*f1a
      i1a0 = (3 - 4*sinhat)/6*fds(y)
      i1b  = coshat/2*(3 - 6*f0b - (2*y + z)*f2b - 4*x*f1b)
      i1b0 = coshat*(4*gds(y) + 1)
      i1cd = - sinhat*x*f1b
      i1ef = (3 - 2*sinhat)/6*x/z*(1 + x/z*dlog(x))
      i2a  = - sinhat/6*x*(1 + 2*f0a + z*f2a) 
     .       - x**2*(3 - 4*sinhat)/12*f1a
      i2b  = (1 - 2*sinhat)/8*x*(3 - 2*f0b - z*f2b - 2*x*f1b)
      i2cd = (3 - 2*sinhat)/12*x*x/z*(1 - z/x + x/z*dlog(x))

      fb = i1a - i1a0 + i1b - i1b0 + i1cd + i1ef + i2a + i2b + i2cd

      dkbtop = - alphat/4/pi1/sinhat*fb

      at = alphas(mt)/pi1
      
C   OLD QCD correction for the subleading mt term, from A. Kwiatkowski 
C   and M. Steinhauser, yellow report 95-03, eq. (10) p. 281, however
C   with rattz2 --> rattw2, following S. Peris and A. Santamaria CERN/TH/95-21.
C
C         dkbtop = dkbtop + alphat/4/pi1/sinhat*
C     .        7/243.d0*(1 + ratzw2/2)*dlog(rattw2)*alphas(mw)/pi1
C  
C   NEW QCD correction for the subleading mt term, from 
C   R. Harlander, T. Seidensticker, and M. Steinhauser, hep-ph/9712228:

      if (fasmt0.eqv..true.) then
         dkbtop = dkbtop + at*alphat/pi1/sinhat*
     .            (0.38d0 + 0.01d0*dlog(x)
     .          + (2.69d0 - 1.50d0*dlog(x)                    )/x
     .          + (0.82d0 - 3.70d0*dlog(x) - 0.09d0*dlog(x)**2)/x**2
     .          - (5.63d0 + 1.44d0*dlog(x) + 0.20d0*dlog(x)**2)/x**3)
      endif

C   If ffermi = true, we substract the leading mt^2 term and add it instead
C   later on, using the Fermi constant G_F.
C   QCD correction for the leading mt**2 term (J. Fleischer, F. Jegerlehner, 
C   P. Raczka and O. V. Tarasov, PLB 293 (1992) 437.)

      xasmt2 = 0.d0
      if (fasmt2.eqv..true.) xasmt2 = (pi2 - 8)/3

      dkbtop = dkbtop - at*alphat/8/pi1/sinhat*rattw2*xasmt2

      if (ffermi.eqv..true.) then
         xt = gf*mt**2/2/pi2/dsqrt(2.d0)
         if (fa2mt4.eqv..false.) dkbtop = dkbtop +
     .        (xt/2 - alphat/8/pi1/sinhat*rattw2)*(1 - at*xasmt2)
      endif

      if (fa2mt4.eqv..true.) then
         dkbtop = dkbtop - alphat/8/pi1/sinhat*rattw2*(1 - at*xasmt2)
         if (ffermi.eqv..false.) xt = alphat/pi1/4/sinhat*rattw2
      endif

C   Leading (alpha mt**2)**2 term with full M_H dependence as comupted by
C   R. Barbieri, M. Beccaria, P. Ciafaloni, G. Curci and A. Vicere, 
C   PLB 288 (1992) 95 and NPB 409 (1993) 105; we use the expressions from
C   J. Fleischer, O.V. Tarasov and F. Jegerlehner, PLB 319 (1993) 249.

      taub2 = 0.d0
      if (fa2mt4.eqv..true.) then
         aa = 1/ratth2
         taub2 = 9 - 13/4.d0*aa - 2*aa**2 - aa*dlog(aa)/4*(19 + 6*aa)
     .         - (aa*dlog(aa))**2/4*(7 - 6*aa) 
     .         - (1 + 14*aa**2 - 12*aa**3)/24*pi2 
     .         + (aa/2 - 2)*dsqrt(aa)*gftj(aa) 
     .         + (aa - 1)**2*(4*aa - 7/4.d0)*fftj0(aa) 
     .         - (aa**3 - 33/4.d0*aa**2 + 18*aa -7)*fftj1(aa)
      endif

      taub = - xt/2*(1 + xt/4*taub2 - at*xasmt2)

      if (mt.le.mz/2) then
         dkbtop = 0.d0
           taub = 0.d0
      endif

      return
      end
      
      
      complex*16 function fds(x)

C   G. Degrassi and A. Sirlin, NPB 352 (1991) 342.

      implicit none
      double precision x,xplus
      complex*16 spence
      include 'ff.h'
      include 'common.f'

      xplus = 1 + x
        fds = 0.5d0 + (3 + 2/x)*(1 - dlog(x)) + (xplus/x)**2
     .        * (2*spence((1,0)/xplus) - 2*zeta2 + dlog(xplus)**2)
     .        + (0,1)*pi1*(2/x + 3 - 2*(xplus/x)**2*dlog(xplus))

      return
      end


      double precision function gds(x)

C   G. Degrassi and A. Sirlin, NPB 352 (1991) 342.

      implicit none
      double precision x,lambda

        gds = 9/8.d0 + 0.5d0/x + (0.5d0 + 1/x)*(lambda(x) - 1) 
     .        - (4 + 2/x)/x*datan(sqrt(x/(4 - x)))**2

      return
      end


      double precision function fftj0(x)

C   J. Fleischer, O.V. Tarasov and F. Jegerlehner, PLB 319 (1993) 249.

      implicit none
      double precision x
      complex*16 spence
      include 'ff.h'
      include 'common.f'

      if (x.lt.1.d0) then
         fftj0 = zeta2 - spence((1,0)*x) - dlog(x)*dlog(1 - x)
      else
         if (x.eq.1.d0) then
            fftj0 = (0.d0,0.d0)
         else
            fftj0 = - zeta2 + spence((1,0)/x) 
     .            - dlog(x)*(dlog(x - 1) - dlog(x)/2)
         endif
      endif

      return
      end


      double precision function fftj1(x)

C   J. Fleischer, O.V. Tarasov and F. Jegerlehner, PLB 319 (1993) 249.

      implicit none
      double precision x,y,zeta,phi
      complex*16 spence
      include 'ff.h'
      include 'common.f'

      y = 4/x

      if (x.gt.4.d0) then
         zeta = (dsqrt(1 - y) - 1)/(dsqrt(1 - y) + 1)
         fftj1 = - (zeta2 + 2*spence((1,0)*zeta) + dlog( - zeta)**2/2)
     .           /dsqrt(1 - y)
      else
         if (x.eq.4.d0) then
            fftj1 = - 4*dlog(2.d0)
         else
            phi = 2*dasin(dsqrt(x/4))
            fftj1 = - 2*dimag(spence(exp((0,1)*phi)))/dsqrt(y - 1)
         endif
      endif

      return
      end


      double precision function gftj(x)

C   J. Fleischer, O.V. Tarasov and F. Jegerlehner, PLB 319 (1993) 249.

      implicit none
      double precision x,y,zeta,phi
      include 'common.f'

      if (x.ge.4.d0) then
         y = 4/x
         zeta = (dsqrt(1 - y) - 1)/(dsqrt(1 - y) + 1)
         gftj = dsqrt(x - 4)*dlog( - zeta)
      else
         phi = 2*dasin(dsqrt(x/4))
         gftj = dsqrt(4 - x)*(pi1 - phi)
      endif

      return
      end
