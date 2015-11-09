      program sintheta

      implicit none
      integer i,j,ssteps
      double precision chi2,lstmin,lstinc,norm,cl,cl1
      double precision lnsint(0:1400),probst(0:1400)
      external fcn
      external chi2
      include '../core/common.f'

      call mintio(5,6,9)

      open (1,file='sinth_zprime.out',status='unknown')
      open (2,file='sinth_plot.out',status='unknown')
      open (8,file='sinth_cl.out',status='unknown')
      open (6,file='/dev/null',status='unknown')

      fwrite = .false.
      fsplot = .false.
      flprob = .true.
      flagmh = .true.
      flagmt = .true.
      flagmc = .true.
      flagal = .true.
      flagS  = .true.
      flagT  = .true.
      flgrho = .true.
      fkappa = .true.
      fsinth = .false.
      fzprim = .true.

      lstmin =  -11.5d0
      lstinc =  1.d0
      ssteps =  7

      norm = 0.d0
      do 10 i = 0, ssteps
         lnsint(i) = lstmin + i*lstinc
         sinth = dexp(lnsint(i))
         open (5,file='sinth.dat',status='old')
         call minuit(fcn,chi2)
         close (5)
         probst(i) = prob
         norm = norm + prob
         write(2,*) lnsint(i)/dlog(10.d0),-2.d0*dlog(probst(i))
 10   continue

      write(8,*)
      write(8,*) 'posterior quantiles:'
      write(8,*)

      cl = 0.d0
      do 20 i = 0, ssteps
         probst(i) = probst(i)/norm
                cl = cl + probst(i)
         write(1,100) i,lnsint(i),dexp(lnsint(i)),1.d2*probst(i),1.d2*cl

         if ((cl.ge.0.05d0).and.(cl1.lt.0.05d0))
     .        write(8,200) ' 5.   % CL: ',
     .             lnsint(i) + lstinc*(cl - 0.05d0)/(cl1 - cl),
     .        dexp(lnsint(i) + lstinc*(cl - 0.05d0)/(cl1 - cl))

         if ((cl.ge.0.15866d0).and.(cl1.lt.0.15866d0)) 
     .        write(8,200) '15.866% CL: ',
     .             lnsint(i) + lstinc*(cl - 0.15866d0)/(cl1 - cl),
     .        dexp(lnsint(i) + lstinc*(cl - 0.15866d0)/(cl1 - cl))

         if ((cl.ge.0.50d0).and.(cl1.lt.0.50d0)) 
     .        write(8,200) '50.   % CL: ',
     .             lnsint(i) + lstinc*(cl - 0.50d0)/(cl1 - cl),
     .        dexp(lnsint(i) + lstinc*(cl - 0.50d0)/(cl1 - cl))

         if ((cl.ge.0.68269d0).and.(cl1.lt.0.68269d0)) 
     .        write(8,200) '68.269% CL: ',
     .             lnsint(i) + lstinc*(cl - 0.68269d0)/(cl1 - cl),
     .        dexp(lnsint(i) + lstinc*(cl - 0.68269d0)/(cl1 - cl))

         if ((cl.ge.0.84134d0).and.(cl1.lt.0.84134d0)) 
     .        write(8,200) '84.134% CL: ',
     .             lnsint(i) + lstinc*(cl - 0.84134d0)/(cl1 - cl),
     .        dexp(lnsint(i) + lstinc*(cl - 0.84134d0)/(cl1 - cl))

         if ((cl.ge.0.90d0).and.(cl1.lt.0.90d0)) 
     .        write(8,200) '90.   % CL: ',
     .             lnsint(i) + lstinc*(cl - 0.90d0)/(cl1 - cl),
     .        dexp(lnsint(i) + lstinc*(cl - 0.90d0)/(cl1 - cl))

         if ((cl.ge.0.95d0).and.(cl1.lt.0.95d0)) 
     .        write(8,200) '95.   % CL: ',
     .             lnsint(i) + lstinc*(cl - 0.95d0)/(cl1 - cl),
     .        dexp(lnsint(i) + lstinc*(cl - 0.95d0)/(cl1 - cl))

         if ((cl.ge.0.99d0).and.(cl1.lt.0.99d0)) 
     .        write(8,200) '99.   % CL: ',
     .             lnsint(i) + lstinc*(cl - 0.99d0)/(cl1 - cl),
     .        dexp(lnsint(i) + lstinc*(cl - 0.99d0)/(cl1 - cl))

         cl1 = cl
 20   continue

 100  format(i4,f9.3,2d12.3,f14.8)
 200  format(a11,f9.3,d12.3)

      stop
      end
