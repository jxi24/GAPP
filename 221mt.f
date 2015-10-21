      program ttomt

      implicit none
      integer i, mtSteps
      double precision mtMin, mtMax, mtDelta
      double precision chisqr, mtProb, norm
      double precision mtMean, mt2Exp

      external fcn
      external chi2
      include 'common.f'

      call mintio(5,6,7)

      flgfitx = .true.
      flgtph  = .true.
      flgs2b  = .true.
      flfout  = .false.
      fwrite  = .false.
      flprob  = .false.
      flagmh  = .true.
      flagmt  = .false.
      flagmc  = .true.
      flagal  = .true.
      flagS   = .true.
      flagT   = .true.
      flgrho  = .true.
      fkappa  = .true.
      fzprim  = .true.
      fsinth  = .true.

      open (6,file='/dev/null',                      status='unknown')
      open (20,file='221chi2/mt_NP/chi2_mt_nu-d.out', status='unknown')
      open (21,file='221chi2/mt_NP/mt_nu-d.out',      status='unknown')

      mtSteps = 200
      mtMin   =  50.0d0
      mtMax   = 250.0d0
      mtDelta = (mtMax-mtMin) / mtSteps

      norm   = 0.0d0
      mtMean = 0.0d0
      mt2Exp = 0.0d0

      do 10 i = 0, mtSteps

         mt = mtMin + i*mtDelta

         open (5,file='221chi2/mt_NP/plot_nu-d.dat',status='old')
         call minuit(fcn,chi2)
         close (5)

         chisqr = prob
         mtProb = dexp( - chisqr/2.0d0)
         norm = norm + mtProb

         mtMean = mtMean + mt    * mtProb
         mt2Exp = mt2Exp + mt**2 * mtProb

C        write(20,100) lnmh, '\t', mh, '\t', chisqr
         write(20,105) mt, '\t', chisqr

 10   continue
   
      mtMean = mtMean / norm
      mt2Exp = mt2Exp / norm

      write(21,*) 
     . 'Exp. Value\tsqrt(Var)\tmt (low)\tmt (central)\tmt (high)\n'

      write(21,110) 
     . mtMean, '\t',
     . sqrt(mt2Exp - mtMean**2), '\t',
     . mtMean - sqrt(mt2Exp - mtMean**2), '\t', 
     . mtMean, '\t',
     . mtMean + sqrt(mt2Exp - mtMean**2)
    
C100  format(f8.4,A,f8.4)
 105  format(f8.4,A,f8.4,A,f8.4)
 110  format(f8.4,A,f8.4,A,f8.4,A,f8.4,A,f8.4)
     
      stop
      end
