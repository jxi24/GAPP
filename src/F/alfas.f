      double precision function alphas(x)

C The used beta function coefficients are from K'uhn et al. 
C (yellow report 95-03, eq. (42) p. 186) and the matching coefficient (11/72)
C is from ibid. eq. (95) p. 196. (Bernreuther and Wetzel in NPB 197 (1982) 228
C seem to give (7/72); the numerical difference is, however, very small.)

      implicit none
      double precision x,arun,cmatch,c3m
      include 'common.f'

      cmatch = 11/72.d0
      if (mu0.lt.mt) then
         if (x.ge.mt) then
            alphas = arun(5,mu0,mt,alfas0)
            if (f4lqcd.eqv..false.) then
               alphas = alphas/(1 + cmatch*(alphas/pi1)**2)
            else
               alphas = alphas/(1 + cmatch*(alphas/pi1)**2
     .                            + c3m(5)*(alphas/pi1)**3)
            endif
            alphas = arun(6,mt,x,alphas)
         else
            if (x.lt.mb) then
               alphas = arun(5,mu0,mb,alfas0)
               if (f4lqcd.eqv..false.) then
                  alphas = alphas*(1 + cmatch*(alphas/pi1)**2)
               else
                  alphas = alphas*(1 + cmatch*(alphas/pi1)**2
     .                               + c3m(4)*(alphas/pi1)**3)
               endif
               if (x.lt.mc) then
                  alphas = arun(4,mb,mc,alphas)
                  if (f4lqcd.eqv..false.) then
                     alphas = alphas*(1 + cmatch*(alphas/pi1)**2)
                  else
                     alphas = alphas*(1 + cmatch*(alphas/pi1)**2
     .                                  + c3m(3)*(alphas/pi1)**3)
                  endif
                  alphas = arun(3,mc,x,alphas)
               else
                  alphas = arun(4,mb,x,alphas)
               endif
            else
               alphas = arun(5,mu0,x,alfas0)
            endif
         endif
      else
         if (x.lt.mt) then
            alphas = arun(6,mu0,mt,alfas0)
            if (f4lqcd.eqv..false.) then
               alphas = alphas*(1 + cmatch*(alphas/pi1)**2)
            else
               alphas = alphas*(1 + cmatch*(alphas/pi1)**2
     .                            + c3m(5)*(alphas/pi1)**3)
            endif
            if (x.lt.mb) then
               alphas = arun(5,mt,mb,alphas)
               if (f4lqcd.eqv..false.) then
                  alphas = alphas*(1 + cmatch*(alphas/pi1)**2)
               else
                  alphas = alphas*(1 + cmatch*(alphas/pi1)**2
     .                               + c3m(4)*(alphas/pi1)**3) 
               endif
               if (x.lt.mc) then
                  alphas = arun(4,mb,mc,alphas)
                  if (f4lqcd.eqv..false.) then
                     alphas = alphas*(1 + cmatch*(alphas/pi1)**2)
                  else
                     alphas = alphas*(1 + cmatch*(alphas/pi1)**2
     .                                  + c3m(3)*(alphas/pi1)**3)   
                  endif
                  alphas = arun(3,mc,x,alphas)
               else
                  alphas = arun(4,mb,x,alphas)
               endif
            else
               alphas = arun(5,mt,x,alphas)
            endif
         else
            alphas = arun(6,mu0,x,alfas0)
         endif
      endif

C For comparison, here is the running alphas formula in terms of 
C Lambda_QCD^(5) = 0.233 according to formula (45) in K"uhn et al of 
C CERN yellow report 95-03, p. 187. 
C
C      l = dlog((x/0.233d0)**2)
C      alphas = pi1/beta0/l*(1 - beta1/beta0**2*dlog(l)/l +
C     .         1/(beta0*l)**2*((beta1/beta0)**2*
C     .         (dlog(l)**2 - dlog(l) - 1) + beta2/beta0))

      return 
      end


      double precision function arun(nf,mu1,mu2,arun0)

      implicit none
      double precision mu1,mu2,arun0,a0,abetal,ab1,ab2
      double precision beta0,b1,b2,b3,fbeta0,fbeta1,fbeta2,fbeta3
      integer nf
      include 'common.f'      
  
      beta0 = fbeta0(nf)      
      a0 = arun0/pi1
      abetal = a0*beta0*dlog((mu2/mu1)**2)

      if (f4lqcd.eqv..false.) then
           b1 = fbeta1(nf)/beta0
           b2 = fbeta2(nf)/beta0
         arun = arun0/(1 + abetal + a0*b1*
     .            dlog(1 + abetal + a0*b1*dlog(1 + abetal))
     .        + a0**2*abetal*(b2 - b1**2)/(1 + abetal))
      else
           b1 = fbeta1(nf)/beta0
           b2 = fbeta2(nf)/beta0
           b3 = fbeta3(nf)/beta0
C  5-loop effects to study uncertainty:
C           b3 = b3 + 579.d0/beta0*a0/(1+abetal + a0*b1*
C     .            dlog(1 + abetal + a0*b1*dlog(1 + abetal))
C     .        + a0**2*abetal*(b2 - b1**2)/(1 + abetal))          
          ab1 = a0*b1
          ab2 = a0**2*(b2 - b1**2)
         arun = arun0/(1 + abetal 
     .        + ab1*dlog(1 + abetal + ab1*dlog(1 + abetal 
     .        + ab1*dlog(1 + abetal)) + abetal*ab2/(1 + abetal))
     .        + ab2*(1 - 1/(1 + abetal + ab1*dlog(1 + abetal)))
     .        + a0**3*(b3/2 - b1*b2 + b1**3/2)*(1 - 1/(1 + abetal)**2))
      endif
           
      return
      end


      double precision function fbeta0(nf)
      integer nf
      fbeta0 = (11 - 2/3.d0*nf)/4
      return 
      end


      double precision function fbeta1(nf)
      integer nf
      fbeta1 = (102 - 38/3.d0*nf)/16
      return 
      end


      double precision function fbeta2(nf)
      integer nf
      fbeta2 = (2857/2.d0 - 5033/18.d0*nf + 325/54.d0*nf**2)/64
      return 
      end


      double precision function fbeta3(nf)
      integer nf
      include 'common.f'
      fbeta3 = (  149753/  6.d0 +       3564*zeta3 
     .       -  (1078361/162.d0 + 6508/27.d0*zeta3)*nf 
     .       +  (  50065/162.d0 + 6472/81.d0*zeta3)*nf**2 
     .       +      1093/729.d0                    *nf**3)/256
      return 
      end


      double precision function c3m(nf)
      integer nf
      include 'common.f'
      c3m = - 82043/27648.d0*zeta3 + 564731/124416.d0 - 2633/31104.d0*nf
      return 
      end
