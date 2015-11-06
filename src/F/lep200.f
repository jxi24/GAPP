      subroutine lep200(xs200,afb200)

C  Calculates cross sections and asymmetries at LEP2 energies. ISR-FSR 
C  interference is not included in correspondence with definition 1 in 
C  LEP2FF/99-01. Formulas are tree level plus leading QED and (massless) 
C  QCD corrections.

      implicit none
      integer i,f
      double precision xs200(2,0:10),afb200(2,0:10),ddilog,hbarc2,intcut
      double precision cum0a2,cum0a3,al,alfaem,dgamma,alphas,as
      double precision ratzs,sqrts,logsme,beta,beta2,fac1,fac2,var
      double precision cr,cq,ci,crf,cqf,cif,arf,aif
      include 'common.f'

      hbarc2 = 0.38937966d9
      intcut = 0.2775d0
      cum0a2 =        85/8.d0 -     23*zeta3/3
      cum0a3 = 372841/2592.d0 -  15961*zeta3/108  + 575*zeta5/18 
     .                        -    529*zeta2/72

      al = alpha/pi1

      do 100 i = 1,2
         if (i.eq.1) sqrts = 182.69d0
         if (i.eq.2) sqrts = 188.63d0
         ratzs = mz2/sqrts/sqrts
         
         call alfahat(sqrts,dgamma,alfaem)
         as = alphas(sqrts)/pi1
         
         logsme = 2*dlog(sqrts/me)
         beta   = 2*al*(logsme - 1)
         beta2  = beta*beta
         fac1   = 4*pi1/3*alfaem**2
         fac2   = intcut**beta*(1 + al*(3*logsme/2 + pi2/3 - 2))
         var    = intcut/(1 - ratzs)
         
         cr = fac1*(fac2/(1 - ratzs)**2*(1 + beta*(beta-1)*dlog(1 - var)
     .      + beta*ratzs*var/(1 - var) - beta2*ddilog(var)) 
     .      - beta/2*(intcut - (1 + 2*ratzs)*dlog(1 - var) 
     .      + ratzs*(1 + ratzs)*var/(1 - intcut - ratzs)))
         ci = fac1*(fac2/(1 - ratzs)*(1 - beta*dlog(1 - var) - beta2*
     .        ddilog(var)) - beta/2*(intcut - (1+ratzs)*dlog(1 - var)))
         cq = fac1*(fac2*(1 - beta*dlog(1 - intcut) 
     .      - beta2*ddilog(intcut)) - beta/2*(intcut - dlog(1- intcut)))

         xs200(i,10) = 0.d0
         do 50 f = 0, 9
            crf =   cr*(v(1)**2 + a(1)**2)*(v(f)**2 + a(f)**2)
            cif = 2*ci*v(1)*v(f)*q(1)*q(f)
            cqf =   cq*q(f)**2
            arf = 4*cr*v(1)*v(f)*a(1)*a(f)
            aif = 2*ci*q(1)*q(f)*a(1)*a(f)

            crf = (v(1)**2 + a(1)**2)*(v(f)**2 + a(f)**2)/(1 - ratzs)**2
     .           *4*pi1/3*alfaem**2
            cif = 2*v(1)*v(f)*q(1)*q(f)/(1 - ratzs)*4*pi1/3*alfaem**2
            cqf =   q(f)**2*4*pi1/3*alfaem**2
            arf = 4*v(1)*v(f)*a(1)*a(f)/(1 - ratzs)**2*4*pi1/3*alfaem**2
            aif = 2*q(1)*q(f)*a(1)*a(f)/(1 - ratzs)*4*pi1/3*alfaem**2
            
            if (f.ne.6) then
                xs200(i,f) = nc(f)*(crf + cif + cqf)/sqrts**2*hbarc2*
     .                       (1 + 3*alfaem/4/pi1*q(f)**2)
               afb200(i,f) = 3*(arf + aif)/(crf + cif + cqf)/4
               if (f.ge.4) then
                  xs200(i,f) = xs200(i,f)*
     .                         (1 + as + as*as*(cum0a2 + as*cum0a3)) 
                  xs200(i,10) = xs200(i,10) + xs200(i,f)
               endif
            else
                xs200(i,f) = 0.d0
               afb200(i,f) = 0.d0
            endif
 50      continue
 100  continue
      
      return
      end




