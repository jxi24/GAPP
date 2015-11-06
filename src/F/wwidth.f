      subroutine wwprod(gammaw)

      implicit none
      double precision gammaw(7),alphas,aw
      include 'common.f'

      aw = alphas(mw)/pi1

C Kai014 ---------------------------------------
C     Corrections to the W width
C ----------------------------------------------
C     Loop factor:
C
C     *(1 - 0.00355d0)  
C ----------------------------------------------

      if( modtype.eq.0 ) then

        gammaw(1) = gf*mw**3/6/dsqrt(2.d0)/pi1
        gammaw(3) = gammaw(1)
        gammaw(4) = gammaw(1)

      endif

      if( modtype.eq.1 ) then
      
         gammaw(1) = (gf*(-3*fitcph**2*kkcc + 3*fits2b*kkcc + 
     -      2*fitx*(-kkcc + kkss))*Sqrt(pi1))/
     -  (24.*2**0.25*fitx*((gf*kkss)/alphat)**1.5*(-kkcc + kkss))

         gammaw(3) = gammaw(1)
         gammaw(4) = gammaw(1)

      endif
           
      if( modtype.eq.2 ) then
      
        gammaw(1) = (gf*(-3*fitcph**2*kkcc + 6*fits2b*kkcc + 
     -      8*fitx*(-kkcc + kkss))*Sqrt(pi1))/
     -  (96.*2**0.25*fitx*((gf*kkss)/alphat)**1.5*(-kkcc + kkss))

        gammaw(3) = gammaw(1)
        gammaw(4) = gammaw(1)

      endif

      if( modtype.eq.3 ) then
         
         gammaw(1) = (alphat**1.5*(4*fitx + 
     -      fitsph**2*(-7 + 3/(1 - 2*kkss)))*dsqrt(pi1))/
     -  (48.*2**0.25*fitx*dsqrt(gf*kkss**3))

         gammaw(3) = gammaw(1)

         gammaw(4) = (alphat**1.5*(fitx*(2 - 4*kkcc) + 
     -      fitsph*(4 - 8*kkcc + fitsph*(-5 + 7*kkcc)))*
     -    dsqrt(pi1))/(24.*2**0.25*fitx*(-1 + kkcc)*(-1 + 2*kkcc)*
     -    dsqrt(gf - gf*kkcc))
 
      endif

      if( modtype.eq.4 ) then

         gammaw(1) = (alphat**1.5*(fitx*(2 - 4*kkcc) + 
     -   fitsph*(4 - 8*kkcc + fitsph*(-5 + 7*kkcc)))*
     -   dsqrt(pi1))/ (24.*2**0.25*fitx*(-1 + kkcc)*(-1 + 2*kkcc)*
     -   dsqrt(gf - gf*kkcc))

         gammaw(3) = (alphat**1.5*(4*fitx + 
     -   fitsph**2*(-7 + 3/(1 - 2*kkss)))*dsqrt(pi1))/
     -   (48.*2**0.25*fitx*dsqrt(gf*kkss**3))

         gammaw(4) = gammaw(1)

      endif

      gammaw(1) = gammaw(1)*(1 - 0.00355d0)
      gammaw(3) = gammaw(3)*(1 - 0.00355d0)
      gammaw(4) = gammaw(4)*(1 - 0.00355d0)

C Kai014 ---------------------------------------
      gammaw(2) = gammaw(1)      
C Kai014 ---------------------------------------
C     gammaw(3) = gammaw(1) - 0.2d-3
      gammaw(3) = gammaw(3) - 0.2d-3
C     gammaw(4) = gammaw(1)*3*(1 + aw + 1.409d0*aw**2 - 12.77d0*aw**3)
      gammaw(4) = gammaw(4)*3*(1 + aw + 1.409d0*aw**2 - 12.77d0*aw**3)
C ----------------------------------------------
      gammaw(5) = gammaw(4) - 0.1d-3
      gammaw(6) = gammaw(4) + gammaw(5)
      gammaw(7) = gammaw(1) + gammaw(2) + gammaw(3) + gammaw(6)

      return
      end
