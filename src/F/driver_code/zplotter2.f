      program zplotter2

C  Computes the minimal Higgs sector constraints for extra Z's.
C  make sure that flagzp and flgzph in chi2.f are both set to .true.
C  and that the lines
C
C         mzp = dexp(xval(16)) 
C       sinth = xval(17)
C
C  in double precision function chi2 are commented out. Also set
C
C        flagzp = zprime
C

      implicit none
      double precision chi2,chimin,delchi,chiref,prob1!,prob
      double precision mzpmin,mzpinc,sinmin,sininc!,sigma
      integer i,j,ssteps,zsteps!,zprime
      external fcn
      external chi2
      include '../core/common.f'

CC--K      call mintio(7,6,9)
      open (6,file='/dev/null',status='unknown')

      fwrite = .false.
      flprob = .false.
      flagmh = .true.
      flagmt = .true.
      flagmc = .true.

      zprime = 3
      sigma = 1.d10

      open (zprime,file='eta3.out',status='new')

      mzpmin = 300.d0
      mzpinc = 20.d0
      zsteps = 110

      do 10 j = 0, zsteps
         mzp = mzpmin + j*mzpinc
         open (7,file='zplotter.dat',status='old')
CC--K         call minuit(fcn,chi2)
         close (7)
         write(zprime,100) sinth,mzp
 10   continue  

      close(zprime)
      
 100  format(f8.6,f8.1)

      stop
      end

