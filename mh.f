      program higgsbound

      implicit none
      integer i,j,hsteps
      double precision chi2,lmhmin,lmhinc,norm,cl,cl1
      double precision meanln,exvln2,meanmh,exvmh2
      double precision lnmh(0:1000),probmh(0:1000)
      external fcn
      external chi2
      include 'common.f'

      call mintio(7,8,9)

      open (1,file='higgs/mh.out',status='unknown')
      open (2,file='higgs/mh_plot.out',status='unknown')
      open (5,file='higgs/mh_cl.out',status='unknown')
      open (6,file='/dev/null',status='unknown')

      fwrite = .false.
      flprob = .true.
      flagmh = .false.
      flagmt = .true.
      flagmc = .true.
      flagal = .true.
      flagS  = .true.
      flagT  = .true.
      flgrho = .true.
      fkappa = .true.
      fsinth = .false.
      fzprim = .false.

      lmhmin = 4.600d0
      lmhinc = 0.002d0
      hsteps = 1000

c      lmhmin = 1.600d0
c      lmhinc = 0.006d0
c      hsteps = 900

C      lmhmin = 1.600d0
C      lmhinc = 0.008d0
C      hsteps = 938

      norm = 0.d0
      do 11 i = 0, hsteps
         lnmh(i) = lmhmin + i*lmhinc
         mh = dexp(lnmh(i))
         open (7,file='higgs/mh.dat',status='old')
         call minuit(fcn,chi2)
         close (7)
         probmh(i) = prob
         norm = norm + prob
         write(2,*) lnmh(i)/dlog(10.d0),-2.d0*dlog(probmh(i))
c        write(2,*) mh,probmh(i)/mh
 11   continue

      write(5,*)
      write(5,*) 'posterior quantiles:'
      write(5,*)

          cl = 0.d0
      meanln = 0.d0
      meanmh = 0.d0
      exvln2 = 0.d0
      exvmh2 = 0.d0
      do 20 i = 0, hsteps
         probmh(i) = probmh(i)/norm
             cl = cl + probmh(i)
         meanln = meanln +      lnmh(i)   *probmh(i)
         meanmh = meanmh + dexp(lnmh(i))  *probmh(i)
         exvln2 = exvln2 +      lnmh(i)**2*probmh(i)
         exvmh2 = exvmh2 + dexp(lnmh(i)*2)*probmh(i)

         write(1,100) i,lnmh(i),dexp(lnmh(i)),1.d2*probmh(i),1.d2*cl

         if ((cl.ge.0.05d0).and.(cl1.lt.0.05d0)) 
     .        write(5,200) ' 5.   % CL: ',
     .             lnmh(i) + lmhinc*(cl - 0.05d0)/(cl1 - cl),
     .        dexp(lnmh(i) + lmhinc*(cl - 0.05d0)/(cl1 - cl))

         if ((cl.ge.0.15866d0).and.(cl1.lt.0.15866d0)) 
     .        write(5,200) '15.866% CL: ',
     .             lnmh(i) + lmhinc*(cl - 0.15866d0)/(cl1 - cl),
     .        dexp(lnmh(i) + lmhinc*(cl - 0.15866d0)/(cl1 - cl))

         if ((cl.ge.0.50d0).and.(cl1.lt.0.50d0)) 
     .        write(5,200) '50.   % CL: ',
     .             lnmh(i) + lmhinc*(cl - 0.50d0)/(cl1 - cl),
     .        dexp(lnmh(i) + lmhinc*(cl - 0.50d0)/(cl1 - cl))

         if ((cl.ge.0.68269d0).and.(cl1.lt.0.68269d0)) 
     .        write(5,200) '68.269% CL: ',
     .             lnmh(i) + lmhinc*(cl - 0.68269d0)/(cl1 - cl),
     .        dexp(lnmh(i) + lmhinc*(cl - 0.68269d0)/(cl1 - cl))

         if ((cl.ge.0.84134d0).and.(cl1.lt.0.84134d0)) 
     .        write(5,200) '84.134% CL: ',
     .             lnmh(i) + lmhinc*(cl - 0.84134d0)/(cl1 - cl),
     .        dexp(lnmh(i) + lmhinc*(cl - 0.84134d0)/(cl1 - cl))

         if ((cl.ge.0.90d0).and.(cl1.lt.0.90d0)) 
     .        write(5,200) '90.   % CL: ',
     .             lnmh(i) + lmhinc*(cl - 0.90d0)/(cl1 - cl),
     .        dexp(lnmh(i) + lmhinc*(cl - 0.90d0)/(cl1 - cl))

         if ((cl.ge.0.95d0).and.(cl1.lt.0.95d0)) 
     .        write(5,200) '95.   % CL: ',
     .             lnmh(i) + lmhinc*(cl - 0.95d0)/(cl1 - cl),
     .        dexp(lnmh(i) + lmhinc*(cl - 0.95d0)/(cl1 - cl))

         if ((cl.ge.0.99d0).and.(cl1.lt.0.99d0)) 
     .        write(5,200) '99.   % CL: ',
     .             lnmh(i) + lmhinc*(cl - 0.99d0)/(cl1 - cl),
     .        dexp(lnmh(i) + lmhinc*(cl - 0.99d0)/(cl1 - cl))

         cl1 = cl
 20   continue

      write(5,*)
      write(5,*) 'E(lnmh,mh) sqrt(var)    mh_min mh_central mh_max'
      write(5,*)

      write(5,300) meanln,sqrt(exvln2 - meanln**2),
     .     dexp(meanln - sqrt(exvln2 - meanln**2)),dexp(meanln),
     .     dexp(meanln + sqrt(exvln2 - meanln**2))

      write(5,300) meanmh,sqrt(exvmh2 - meanmh**2),
     .     meanmh - sqrt(exvmh2 - meanmh**2),meanmh,
     .     meanmh + sqrt(exvmh2 - meanmh**2)

 100  format(i4,f8.3,f10.2,d12.5,f14.8)
 200  format(a11,2f9.4)
 300  format(1x,f9.4,f10.4,'   ',3f9.4)

      stop
      end
