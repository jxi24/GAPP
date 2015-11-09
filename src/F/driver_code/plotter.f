      program mwmtplot

      implicit none
      integer i,j,tsteps,hsteps
      double precision chi2,chiref,chimin,delchi,prob1,mhmin,mhinc
      double precision mtval,mterr,mwval,mwerr,mtmin,mtinc,mw1
      external fcn
      external chi2
      include '../core/common.f'

      open ( 8,file='mwmtplot/mwmt_117.out',   status='unknown')
      open ( 9,file='mwmtplot/mwmt_200.out',   status='unknown')
      open (10,file='mwmtplot/mwmt_300.out',   status='unknown')
      open (11,file='mwmtplot/mwmt_500.out',   status='unknown')
      open (12,file='mwmtplot/mwmt_direct.out',status='unknown')
      open (13,file='mhmtplot/mhmt_mt.out',    status='unknown')

CC--K      call mintio(5,6,7)

      open (6,file='/dev/null',status='unknown')
      
      fwrite = .false.
      flprob = .false.
      flagmh = .false.
      flagmt = .false.
      flagmc = .true.
      flagal = .false.
      flagS  = .true.
      flagT  = .true.
      flgrho = .true.
      fkappa = .true.
      fzprim = .true.
      fsinth = .true.

        mtmin = 140.00d0
        mtinc =   0.50d0
       tsteps = 100
      
C  enter central values and errors of direct measurements here:

      mtval = 170.930d0
      mterr =   1.930d0
      mwval =  80.398d0
      mwerr =   0.025d0

      dahad3 = 0.005571d0
      do 20 i = 8, 11
         if (i.eq. 8) mh =  117.d0
         if (i.eq. 9) mh =  200.d0
         if (i.eq.10) mh =  300.d0
         if (i.eq.11) mh =  500.d0
         do 10 j = 0, tsteps
            mt = mtmin + j*mtinc
            open (5,file='mwmtplot/plot.dat',status='old')
CC--K            call minuit(fcn,chi2)
            close(5)
            write(i,100) mtp, mw
 10      continue
 20   continue

      dahad3 = 0.005811d0
      do 40 i = 8, 11
         if (i.eq. 8) mh =  117.d0
         if (i.eq. 9) mh =  200.d0
         if (i.eq.10) mh =  300.d0
         if (i.eq.11) mh =  500.d0
         do 30 j = 0, tsteps
            mt = mtmin + (tsteps - j)*mtinc
            open (5,file='mwmtplot/plot.dat',status='old')
CC--K            call minuit(fcn,chi2)
            close(5)
            write(i,100) mtp, mw
 30      continue
 40   continue

      do 50 i = 1, 90
         write(12,100) mtval - sin(i*pi1/45.d0)*mterr, 
     .                 mwval - cos(i*pi1/45.d0)*mwerr
 50   continue

      write(13,100) mtval - mterr, 1.5d0
      write(13,100) mtval - mterr, 7.0d0
      write(13,100) mtval + mterr, 7.0d0
      write(13,100) mtval + mterr, 1.5d0

 100  format(f8.4,f8.4)
      
      stop
      end
