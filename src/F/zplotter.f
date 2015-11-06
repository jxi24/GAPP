      program zplotter

      implicit none
      double precision chi2,chimin,dchi1s,dchi90,cref1s,cref90,prob1
      double precision mzpmin,mzpinc,sinmin,sininc
      double precision eplus(20),eminus(20),eparab(20),globcc(20)
      integer i,j,ssteps,zsteps
      external fcn
      external chi2
      include 'common.f'

      ssteps =  4
      zsteps =  9

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
      fzprim = .true.

      call mintio(5,6,7)
      open (5,file='smfit.dat',status='old')
      open (6,file='/dev/null',status='unknown')
      call minuit(fcn,chi2)
      close(5)

      fsinth = .false.
      fzprim = .false.

      open (1,file='chi_1.out',status='unknown')
      open (2,file='chi_2.out',status='unknown')
      open (8,file='chi_3.out',status='unknown')
      open (9,file='chi_4.out',status='unknown')

      dchi1s =  1.00d0
c      dchi90 =  3.04d0
c      dchi90 =  4.6052d0
      cref1s = prob + dchi1s
      cref90 = prob + dchi90

c      call mnerrs(16,eplus(16),eminus(16),eparab(16),globcc(16))
c      call mnerrs(17,eplus(17),eminus(17),eparab(17),globcc(17))
c      mzpmin = dlog(mzp) + eminus(16)
c      sinmin =     sinth + eminus(17)
c      mzpinc = (eplus(16) - eminus(16))/(2*zsteps - 2)
c      sininc = (eplus(17) - eminus(17))/(2*ssteps - 2)

C  if that fails:

      mzpmin =  6.2146d0
      mzpinc =  0.2d0
      sinmin = -0.0010d0
      sininc =  0.0002d0

      do 40 j = 0, ssteps
         sinth = sinmin + (2*j - 1)*sininc
         prob1 = 1.d4
         do 30 i = 0, zsteps
            mzp = dexp(mzpmin + (2*i - 1)*mzpinc)
            open (5,file='zplotter.dat',status='old')
            call minuit(fcn,chi2)
            close (5)
            if ((prob1.ge.cref1s).and.(prob.le.cref1s))
     .           write(1,100) sinth, dexp(mzpmin + 2*mzpinc*
     .           (i - 1/2.d0 + (prob - cref1s)/(prob1 - prob)))
            if ((prob1.le.cref1s).and.(prob.ge.cref1s))
     .           write(2,100) sinth, dexp(mzpmin + 2*mzpinc*
     .           (i - 1/2.d0 + (prob - cref1s)/(prob1 - prob)))
            if ((prob1.ge.cref90).and.(prob.le.cref90))
     .           write(8,100) sinth, dexp(mzpmin + 2*mzpinc*
     .           (i - 1/2.d0 + (prob - cref90)/(prob1 - prob)))
            if ((prob1.le.cref90).and.(prob.ge.cref90))
     .           write(9,100) sinth, dexp(mzpmin + 2*mzpinc*
     .           (i - 1/2.d0 + (prob - cref90)/(prob1 - prob)))
            prob1 = prob
 30      continue
 40   continue
      
c      prob1 = 1.d4
c      do 25 i = 0, ssteps
c         sinth = sinmin + i*sininc
c         mzp = 2.5d3
c         open (7,file='zplotter.dat',status='old')
c         call minuit(fcn,chi2)
c         close (7)
c         if (((prob1.ge.chiref).and.(prob.le.chiref)).or.
c     .        ((prob1.le.chiref).and.(prob.ge.chiref))) then
c            sinth = sinth + sininc*(prob - chiref)/(prob1 - prob)
c            write(1,100) sinth,mzp
c         endif
c         prob1 = prob
c 25   continue
      
 100  format(f8.6,f8.1)

      stop
      end
