      program ttomh

      implicit none
      integer i, lnSteps
      double precision lnmh, lnmh1, lnMin, lnMax, lnDelta
      double precision chisqr, mhProb, norm
      double precision mhMean, mh2Exp, lnMean, ln2Exp
      double precision chiref, delchi, prob1

      external fcn
      external chi2
      include '../core/common.f'

CC--K      call mintio(5,6,7)

      flgfitx = .true.
      flgtph  = .true.
      flgs2b  = .true.
      flfout  = .false.
      fwrite  = .false.
      flprob  = .false.
      flagmh  = .true.
      flagmt  = .true.
      flagmc  = .true.
      flagal  = .true.
      flagS   = .true.
      flagT   = .true.
      flgrho  = .true.
      fkappa  = .true.
      fzprim  = .true.
      fsinth  = .true.

      open (5,file='smfit.dat',                      status='old')
      open (6,file='/dev/null',                      status='unknown')
      open (20,file='221chi2/mh_SM/chi2_uu-d.out',status='unknown')
      open (21,file='221chi2/mh_SM/mh_uu-d.out',status='unknown')

CC--K      call minuit(fcn,chi2)

      close(5)

      delchi = 3.84146d0 
      chiref = prob + delchi

      flagmh  = .false.

      lnSteps = 10000
      lnMin   = dlog(25.0d0)
      lnMax   = dlog(500.0d0) 
      lnDelta = (lnMax-lnMin) / lnSteps

      norm   = 0.0d0
      mhMean = 0.0d0
      mh2Exp = 0.0d0
      lnMean = 0.0d0
      ln2Exp = 0.0d0

      prob1 = 1.d4

      do 10 i = 0, lnSteps

         lnmh = lnMin + i*lnDelta
         mh = dexp(lnmh)

         open (5,file='221chi2/mh_SM/plot_uu-d.dat',status='old')
CC--K         call minuit(fcn,chi2)
         close (5)

         chisqr = prob
         mhProb = dexp( - chisqr/2.0d0)
         norm = norm + mhProb

         mhMean = mhMean + mh      * mhProb
         lnMean = lnMean + lnmh    * mhProb
         mh2Exp = mh2Exp + mh**2   * mhProb
         ln2Exp = ln2Exp + lnmh**2 * mhProb

C        write(20,100) lnmh, '\t', mh, '\t', chisqr
         write(20,105) mh, '\t', chisqr

         if ((prob1.ge.chiref).and.(prob.le.chiref)) then

         write(21,'(A,f8.4,A,f8.4)') 'Min. MH:\t',
     .   dexp(lnmh + (lnmh - lnmh1)*(prob - chiref)/(prob1 - prob)),
     .   '\t', lnmh + (lnmh - lnmh1)*(prob - chiref)/(prob1 - prob)

         endif

         if ((prob1.le.chiref).and.(prob.ge.chiref)) then

         write(21,'(A,f8.4,A,f8.4)') 'Max. MH:\t',
     .   dexp(lnmh + (lnmh - lnmh1)*(prob - chiref)/(prob1 - prob)),
     .   '\t', lnmh + (lnmh - lnmh1)*(prob - chiref)/(prob1 - prob)

         endif

         prob1 = prob
         lnmh1 = lnmh

 10   continue
   
      mhMean = mhMean / norm
      lnMean = lnMean / norm
      mh2Exp = mh2Exp / norm
      ln2Exp = ln2Exp / norm


      write(21,*) 
     . '\n\nExp. Value\tsqrt(Var)\tMH (low)\tMH (central)\tMH (high)\n'

      write(21,110) 
     . lnMean, '\t',
     . sqrt(ln2Exp - lnMean**2), '\t',
     . dexp(lnMean - sqrt(ln2Exp - lnMean**2)), '\t', 
     . dexp(lnMean), '\t',
     . dexp(lnMean + sqrt(ln2Exp - lnMean**2)), 
     . '\tln(MH)'

      write(21,110) 
     . mhMean, '\t',
     . sqrt(mh2Exp - mhMean**2), '\t',
     . mhMean - sqrt(mh2Exp - mhMean**2), '\t',
     . mhMean, '\t',
     . mhMean + sqrt(mh2Exp - mhMean**2),
     . '\tMH'
    
C100  format(f8.4,A,f8.4,A,f8.4)
 105  format(f10.6,A,f10.6)
 110  format(f8.4,A,f8.4,A,f8.4,A,f8.4,A,f8.4,A)
     
      stop
      end
