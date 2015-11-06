      program ttoplot

      implicit none
      integer i,j,k,xsteps,ysteps,zsteps, fcalls
      double precision chi2,chiref,chigrd,chimin,delchi,tolrnc
      double precision prob1,lnx, lnx1, cph
      double precision lnxmin,lnxmax,lnxinc 
      double precision tphmin,tphmax,tphinc
      double precision cphmin,cphmax,cphinc
      double precision s2bmin, s2bmax, s2binc

      logical scnlnx, scntph, scns2b

      external fcn
      external chi2
      include 'common.f'

CC--K      call mintio(5,6,7)

      flgfitx = .true.
      flgtph  = .true.
      flgs2b  = .true.

      flfout = .false.

      fwrite = .false.
      flprob = .false.

C --- SET HIGGS MASS ------ 05/31/09 -- set Higgs mass by hand.
C
C     flagmh = .false.
C     mh = 250.d0
      flagmh = .true.
C -------------------------------------------------------------

      flagmt = .true.
      flagmc = .true.
      flagal = .true.
      flagS  = .true.
      flagT  = .true.
      flgrho = .true.
      fkappa = .true.
      fzprim = .true.
      fsinth = .true.

      open (5,file='smfit.dat',                       status='old')
      open (6,file='/dev/null',                       status='unknown')
      open (9,file='221plots/grid_fp-d_sm.dat', status='unknown')
      open (20,file='221plots/plot_fp-d_sm.out',status='unknown')
C     open (21,file='221plots/plot_fp-d_2_sm.out', status='unknown')

CC--K      call minuit(fcn,chi2)

      close(5)

      write(20,'(A)') 
     - 'fits2b: \tfittph: \tfitx: \t\tchi2 (old): \tchi2:'
C     write(21,'(A)') 
C    - 'fits2b: \tfittph: \tfitx: \t\tchi2 (old): \tchi2:'
      write(9,'(A,F10.4)') 'Min. chi2: ', prob

C Delta Chi2 for a fit of 3 parameters at 95 % CL.
      delchi = 7.81473d0
C Delta Chi2 for a fit of 2 parameters at 95 % CL.
C     delchi = 5.99146d0

      tolrnc = 8.d1

C --- Use SM minimum chi2  ------ 06/01/09
C     prob = 41.94848
C ----------------------------------------

      chiref = prob + delchi
      chigrd = chiref + tolrnc

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C PRE-RUN TO FIND THE GRID WE WANT TO LOOP OVER CCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      scnlnx = .false.
C     scnlnx = .true.   -- 06/01/09
      scntph = .false.
      scns2b = .false.

      xsteps = 200
      ysteps = 200
      zsteps = 200

      lnxmin =  0.0d0
      lnxmax =  9.5d0

      tphmin =  0.0d0
      tphmax =  1.0d3
C     tphmax =  5.0d0   -- 06/01/09
C     tphmax =  3.0d1

C -- SCAN OVER COS(PHI) -- 06/02/09 
C----------------------------------
      cphmin =  0.0d0
      cphmax =  1.0d0
C ---------------------------------

      s2bmin  =  0.0d0
      s2bmax  =  1.0d0
C     s2bmin  =  bests2b - 0.15d0
C     s2bmax  =  bests2b + 0.15d0
C
C     if (s2bmin.le.0.d0) then 
C
C        s2bmin = 0.0d0
C        s2bmax = 0.3d0
C
C     endif
C
C     if (s2bmax.ge.1.d0) then
C
C        s2bmin = 0.7d0
C        s2bmax = 1.0d0
C
C     endif

      lnxinc  = (lnxmax-lnxmin) / xsteps
      s2binc  = (s2bmax-s2bmin) / ysteps
      tphinc  = (tphmax-tphmin) / zsteps

C FIND BOUNDARIES FOR fitx CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if (scnlnx.eqv..true.) then

         flgfitx = .false.
         lnx = lnxmin
         prob1  = 1.0d10

         do 10 i = 0, xsteps

            fitx = dexp(lnx)               

            open (5,file='221plots/plot_fp-d_sm.dat')

CC--K            call minuit(fcn,chi2)   
            close(5)

C           write(9,'(F10.4,A,F10.4)') fitx, ': ', prob

            if ((prob1.ge.chigrd).and.(prob.le.chigrd)) then

C              lnxmax = dlog(dexp(lnx) + 3000.0d0)
               lnxmin = lnxmin
C              lnxmin = lnx
               goto 11

            endif

C           if ((prob1.le.chiref).and.(prob.ge.chiref)) then
C
C              lnxmax = lnx
C
C           endif

            lnx = lnx + lnxinc
            prob1 = prob

 10      continue
      endif
11

C FIND BOUNDARIES FOR fittph CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if (scntph.eqv..true.) then

C        flgfitx = .true.
         fitx = dexp(lnxmax)

         flgtph = .false.
         fittph = tphmin
         prob1  = 1.0d10

         do 15 i = 0, ysteps

            open (5,file='221plots/plot_fp-d_sm.dat')
CC--K            call minuit(fcn,chi2)   
            close(5)

C           write(9,'(A,F10.4,A,F10.4)') '\t', fittph, ': ', prob

            if ((prob1.ge.chigrd).and.(prob.le.chigrd)) then

               tphmin = fittph

            endif

            if ((prob1.le.chiref).and.(prob.ge.chiref)) then

                tphmax = fittph

            endif

            fittph = fittph + tphinc
            prob1 = prob

 15      continue

      endif

C FIND BOUNDARIES FOR fits2b CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if(scns2b.eqv..true.) then

         flgtph = .true.
         flgs2b = .false.
         fits2b = s2bmin
         prob1  = 1.0d10

         do 20 i = 0, zsteps

            open (5,file='221plots/plot_fp-d_sm.dat')
CC--K            call minuit(fcn,chi2)   
            close(5)

C           write(9,'(A,F10.4,A,F10.4)') '\t\t', fits2b, ': ', prob

            if ((prob1.ge.chigrd).and.(prob.le.chigrd)) then

               s2bmin = fits2b

            endif

            if ((prob1.le.chiref).and.(prob.ge.chiref)) then

               s2bmax = fits2b

            endif

            fits2b = fits2b + s2binc               
            prob1 = prob

 20      continue

      endif

      write(9,'(A,F10.4,A,F10.4)') 
     -   '\n\nfitxmin: ', dexp(lnxmin), '\tfitxmax: ', dexp(lnxmax)
      write(9,'(A,F10.4,A,F10.4)') 
     -       'cphmin:  ', cphmin,       '\tcphmax:  ', cphmax
      write(9,'(A,F10.4,A,F10.4)') 
     -       's2bmin:  ', s2bmin,       '\ts2bmax:  ', s2bmax

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C LOOP OVER THE CONSTRUCTED GRID CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      flgfitx = .false.
      flgtph  = .false.
      flgs2b  = .false.

      fcalls = 0

C -- Original grid settings
      xsteps = 100
      ysteps =  60
      zsteps = 100

C -- Test grid
C     xsteps = 10
C     ysteps =  6
C     zsteps = 10

      lnxinc  = (lnxmax-lnxmin) / xsteps
      s2binc  = (s2bmax-s2bmin) / ysteps
      tphinc  = (tphmax-tphmin) / zsteps
      cphinc  = (cphmax-cphmin) / zsteps

C     write(9,'(A)') 
C    - '\n\nfits2b: \tfittph: \tfitx: \t\tchi2 (old): \tchi2:'

C     do 50 k = 0, zsteps
      do 50 k = 0, zsteps-1

C -- SCAN OVER COS(PHI) -- 06/02/09 
C----------------------------------
C        fittph = tphmin + k*tphinc
         cph = cphmax - k*cphinc
         fittph = (1-cph*cph)/(cph*cph)

         do 40 j = 0, ysteps

C        fits2b = bests2b --  plot only best fit value for fits2b. 05/30/09
         fits2b = s2bmin + j*s2binc

            prob1 = 1.d4

            do 30 i = 0, xsteps

               lnx = lnxmin + i*lnxinc

               fitx = dexp(lnx)

               open (5,file='221plots/plot_fp-d_sm.dat')
CC--K               call minuit(fcn,chi2)
               close(5)

               fcalls = fcalls + 1

             if ((prob1.ge.chiref).and.(prob.le.chiref)) then
             write(20,100) fits2b, '\t', fittph, '\t', 
     -       dexp(lnx + (lnx - lnx1)*(prob - chiref)/(prob1 - prob)),
     -       '\t', prob1, '\t', prob
               exit
               endif

C              if ((prob1.le.chiref).and.(prob.ge.chiref)) then
C
C              write(21,100) fits2b, '\t', fittph, '\t', 
C    -         dexp(lnx + (lnx - lnx1)*(prob - chiref)/(prob1 - prob)),
C    -         '\t', prob1, '\t', prob
C
C              endif
C
C              write(9,100) fits2b, '\t', fittph, '\t', fitx,
C    -                      '\t', prob1, '\t', prob

               prob1 = prob
               lnx1 = lnx
 30         continue
 40      continue
 50   continue

      write(9,'(A,I7,A,I7,A)') '\n\n\tFunction calls:\t', 
     - fcalls, ' (', (xsteps+1)*(ysteps+1)*zsteps, ')'

 100  format(f8.4,A,f12.4,A,f9.4,A,f12.4,A,f12.4) 

      stop
      end
