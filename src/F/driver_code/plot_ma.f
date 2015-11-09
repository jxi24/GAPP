      program maplotter

      implicit none
      double precision chi2,chimin,dchi1s,dchi95,cref95,prob1
      double precision mzpmin,mzpinc,sinmin,sininc
      double precision eplus(20),eminus(20),eparab(20),globcc(20)
      integer i,j,ssteps,zsteps,ima
      external fcn
      external chi2
      common   /ma/ ima,ssteps
      include '../core/common.f'
      
      ssteps = 60
      zsteps = 151

      fwrite = .false.
      fsplot = .false.
      flprob = .false.
      flagmh = .true.
      flagmt = .true.
      flagmc = .true.
      flagal = .true.
      flagS  = .true.
      flagT  = .true.
      flgrho = .true.
      fkappa = .true.
      fsinth = .true.

CC--K      call mintio(5,6,7)

      sinth = 0.d0

      mzpmin =  6.9078d0
      mzpinc =  0.01d0
      dchi95 = 2.7055d0

      open (1,file='ma_1.out',status='unknown')
      open (2,file='ma_2.out',status='unknown')

      do 40 ima = 0, ssteps
         fzprim = .true.
         open (5,file='smfit.dat',status='old')
         open (6,file='/dev/null',status='unknown')
         open (7,file='smfit.out',status='unknown')
CC--K         call minuit(fcn,chi2)
         close(5)
         
         fzprim = .false.      

         cref95 = prob + dchi95

         prob1 = 1.d4
         do 30 i = 0, zsteps
            mzp = dexp(mzpmin + (2*i - 1)*mzpinc)
            open (5,file='zplotter.dat',status='old')
CC--K            call minuit(fcn,chi2)
            close (5)
            if (ima.le.ssteps/2) then
               if ((prob1.ge.cref95).and.(prob.le.cref95))
     .              write(1,100) datan(eps2_L(1)/eps2_L(4)), 
     .              1/1000.d0*dexp(mzpmin + 2*mzpinc*
     .              (i - 1/2.d0 + (prob - cref95)/(prob1 - prob)))
            else
               if ((prob1.ge.cref95).and.(prob.le.cref95))
     .              write(1,100) datan(eps2_L(1)/eps2_L(4)) + pi1, 
     .              1/1000.d0*dexp(mzpmin + 2*mzpinc*
     .              (i - 1/2.d0 + (prob - cref95)/(prob1 - prob)))
            endif
            if ((prob1.le.cref95).and.(prob.ge.cref95))
     .           write(2,100) eps2_L(1), dexp(mzpmin + 2*mzpinc*
     .           (i - 1/2.d0 + (prob - cref95)/(prob1 - prob)))
            prob1 = prob
 30      continue
 40   continue
      
 100  format(f8.5,f8.4)

      stop
      end
