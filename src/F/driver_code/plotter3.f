      program plot

      implicit none
      integer i,j,flplot,xsteps,ysteps
      double precision chi2,chiref,delchi,prob1,mw1,mtp1
      double precision mtmin,mhmin,mtinc,mhinc,Smin,Sinc,Tmin,Tinc
      double precision rhomin,kapmin,rhoinc,kapinc
      double precision eplus(30),eminus(30),eparab(30),globcc(30)
      external fcn
      external chi2
      include '../core/common.f'

C  choose flplot =  1 for M_W vs m_t and       all data
C  choose flplot =  2 for M_W vs m_t and  indirect data
C  choose flplot =  3 for M_H vs m_t and       all data
C  choose flplot =  4 for M_H vs m_t and X-section data
C  choose flplot =  5 for M_H vs m_t and asymmetry data
C  choose flplot =  6 for M_H vs m_t and       DIS data
C  choose flplot =  7 for M_H vs m_t and       M_W data
C  choose flplot =  8 for  T  vs  S  and       all data (M_H =  117 GeV)
C  choose flplot =  9 for  T  vs  S  and       all data (M_H =  340 GeV)
C  choose flplot = 10 for  T  vs  S  and       all data (M_H = 1000 GeV)
C  choose flplot = 11 for  T  vs  S  and X-section data
C  choose flplot = 12 for  T  vs  S  and asymmetry data
C  choose flplot = 13 for  T  vs  S  and       DIS data
C  choose flplot = 14 for  T  vs  S  and       M_W data
C  choose flplot = 15 for  T  vs  S  and       APV data
C  choose flplot = 16 for rho_b vs kappa_b
C
C  choose delchi = errordef in smfit.dat

      flplot = 7
      xsteps = 100
      ysteps = 60

      fwrite = .false.
      flprob = .false.
      flagmh = .true.
      flagmt = .true.
      flagmc = .true.
      flagal = .true.
      flagS  = .true.
      flagT  = .true.
      flgrho = .true.
      fkappa = .true.
      fzprim = .true.
      fsinth = .true.

      call mintio(5,6,7)
      open (5,file='smfit.dat',status='old')
      open (6,file='/dev/null',status='unknown')
      call minuit(fcn,chi2)
      close(5)

      if (flplot.le.7) then
         flagmh = .false.
         flagmt = .false.
      else if (flplot.eq.16) then
         fkappa = .false.
         flgrho = .false.
      else
         flagS  = .false.
         flagT  = .false.
      endif
         
      delchi = 1.d0     ! default
      if (flplot.eq. 1) then
         open (1,file='mwmtplot/mwmt_all.out', status='unknown')
         open (2,file='mwmtplot/mwmt_all2.out',status='unknown')
         delchi = 4.6052d0
      else if (flplot.eq. 2) then
         open (1,file='mwmtplot/mwmt_indirect.out', status='unknown')
         open (2,file='mwmtplot/mwmt_indirect2.out',status='unknown')
      else if (flplot.eq. 3) then
         open (1,file='mhmtplot/mhmt_all.out', status='unknown')
         open (2,file='mhmtplot/mhmt_all2.out',status='unknown')
         delchi = 4.6052d0
      else if (flplot.eq. 4) then
         open (1,file='mhmtplot/mhmt_Xsec.out', status='unknown')
         open (2,file='mhmtplot/mhmt_Xsec2.out',status='unknown')
      else if (flplot.eq. 5) then
         open (1,file='mhmtplot/mhmt_asym.out', status='unknown')
         open (2,file='mhmtplot/mhmt_asym2.out',status='unknown')
      else if (flplot.eq. 6) then
         open (1,file='mhmtplot/mhmt_DIS.out', status='unknown')
         open (2,file='mhmtplot/mhmt_DIS2.out', status='unknown')
      else if (flplot.eq. 7) then
         open (1,file='mhmtplot/mhmt_mw.out', status='unknown')
         open (2,file='mhmtplot/mhmt_mw2.out',status='unknown')
      else if (flplot.eq. 8) then
         open (1,file='STplot/ST_all_117.out', status='unknown')
         open (2,file='STplot/ST_all2_117.out',status='unknown')
         delchi = 4.6052d0
      else if (flplot.eq. 9) then
         open (1,file='STplot/ST_all_340.out', status='unknown')
         open (2,file='STplot/ST_all2_340.out',status='unknown')
         delchi = 4.6052d0
      else if (flplot.eq.10) then
         open (1,file='STplot/ST_all_1000.out', status='unknown')
         open (2,file='STplot/ST_all2_1000.out',status='unknown')
         delchi = 4.6052d0
      else if (flplot.eq.11) then
         open (1,file='STplot/ST_Xsec.out', status='unknown')
         open (2,file='STplot/ST_Xsec2.out',status='unknown')
      else if (flplot.eq.12) then
         open (1,file='STplot/ST_asym.out', status='unknown')
         open (2,file='STplot/ST_asym2.out',status='unknown')
      else if (flplot.eq.13) then
         open (1,file='STplot/ST_DIS.out', status='unknown')
         open (2,file='STplot/ST_DIS2.out',status='unknown')
      else if (flplot.eq.14) then
         open (1,file='STplot/ST_mw.out', status='unknown')
         open (2,file='STplot/ST_mw2.out',status='unknown')
      else if (flplot.eq.15) then
         open (1,file='STplot/ST_qw.out', status='unknown')
         open (2,file='STplot/ST_qw2.out',status='unknown')
      else
         open (1,file='Zbbplot/Zbb.out', status='unknown')
         open (2,file='Zbbplot/Zbb2.out',status='unknown')
         delchi = 4.6052d0
      endif

      chiref = prob + delchi

      if (flplot.le.7) then
         if (flplot.ge.4) then
            mtmin = 110.0d0
            mhmin =   2.0d0
            mtinc =  80.0d0/(9*xsteps + 2)
            mhinc =   5.0d0/(2*ysteps - 2)
         else
            call mnerrs(2,eplus(2),eminus(2),eparab(2),globcc(2))
            call mnerrs(7,eplus(7),eminus(7),eparab(7),globcc(7))
            mtmin = mt + eminus(2)
            mhmin = dlog(mh) + eminus(7)
            mtinc = (eplus(2) - eminus(2))/(9*xsteps + 2)
            mhinc = (eplus(7) - eminus(7))/(2*ysteps - 2)
         endif
         do 40 i = 0, ysteps
            mh = dexp(mhmin + (2*i - 1)*mhinc)
            prob1 = 1.d4
            do 30 j = 0, xsteps
               mt = mtmin + (9*j + 1)*mtinc
               if (flplot.eq.1) open (5,file='mwmtplot/plot1.dat')
               if (flplot.eq.2) open (5,file='mwmtplot/plot2.dat')
               if (flplot.eq.3) open (5,file='mhmtplot/plot3.dat')
               if (flplot.eq.4) open (5,file='mhmtplot/plot4.dat')
               if (flplot.eq.5) open (5,file='mhmtplot/plot5.dat')
               if (flplot.eq.6) open (5,file='mhmtplot/plot6.dat')
               if (flplot.eq.7) open (5,file='mhmtplot/plot7.dat')
               call minuit(fcn,chi2)
               close(5)
               if ((prob1.ge.chiref).and.(prob.le.chiref)) then
                  if (flplot.le.2) write(1,100) mtp,
     .                 mw + (mw - mw1)*(prob - chiref)/(prob1 - prob)
                  if (flplot.ge.3) write(1,100) mtp + 
     .                 (mtp - mtp1)*(prob - chiref)/(prob1 - prob),
     .                 dlog(mh),prob,prob1
               endif
               if ((prob1.le.chiref).and.(prob.ge.chiref)) then
                  if (flplot.le.2) write(2,100) mtp,
     .                 mw + (mw - mw1)*(prob - chiref)/(prob1 - prob)
                  if (flplot.ge.3) write(2,100) mtp +
     .                 (mtp - mtp1)*(prob - chiref)/(prob1 - prob),
     .                 dlog(mh),prob,prob1
               endif
               mtp1 = mtp
               prob1 = prob
               mw1 = mw
 30         continue
 40      continue
      else if (flplot.eq.16) then
         call mnerrs(11,eplus(11),eminus(11),eparab(11),globcc(11))
         call mnerrs(12,eplus(12),eminus(12),eparab(12),globcc(12))
         rhomin = Brho   + eminus(11)
         kapmin = Bkappa + eminus(12)
         rhoinc = (eplus(11) - eminus(11))/(2*ysteps - 2)
         kapinc = (eplus(12) - eminus(12))/(9*xsteps + 2)
         do 46 j = 0, xsteps
            Bkappa = kapmin + (9*j + 1)*kapinc
            prob1 = 1.d4
            do 45 i = 0, ysteps
               Brho = rhomin + (2*i - 1)*rhoinc
               open (5,file='STplot/plot16.dat')
               call minuit(fcn,chi2)
               close(5)
               if ((prob1.ge.chiref).and.(prob.le.chiref)) write(1,100) 
     .              rhomin + 2*rhoinc*(i - 1/2.d0 
     .                     + (prob - chiref)/(prob1 - prob)), Bkappa
               if ((prob1.le.chiref).and.(prob.ge.chiref)) write(2,100) 
     .              rhomin + 2*rhoinc*(i - 1/2.d0 
     .                     + (prob - chiref)/(prob1 - prob)), Bkappa
               prob1 = prob
 45         continue
 46      continue
      else
         if ((flplot.eq.12).or.(flplot.ge.14)) then
            Tmin = -1.55d0
            Smin = -1.60d0
            Tinc =  3.10d0/(9*xsteps + 2)
            Sinc =  3.20d0/(2*ysteps - 2)
         else
            call mnerrs(8,eplus(8),eminus(8),eparab(8),globcc(8))
            call mnerrs(9,eplus(9),eminus(9),eparab(9),globcc(9))
            Tmin = Tpar + eminus(8)
            Smin = Spar + eminus(9)
            Tinc = (eplus(8) - eminus(8))/(9*xsteps + 2)
            Sinc = (eplus(9) - eminus(9))/(2*ysteps - 2)
         endif
         do 60 j = 0, xsteps
            Tpar = Tmin + (9*j + 1)*Tinc
C            Tpar = 18.d0/452.d0*3.1d0
            prob1 = 1.d4
            do 50 i = 0, ysteps
               Spar = Smin + (2*i - 1)*Sinc
               if (flplot.eq. 8) open (5,file='STplot/plot8.dat')
               if (flplot.eq. 9) open (5,file='STplot/plot9.dat')
               if (flplot.eq.10) open (5,file='STplot/plot10.dat')
               if (flplot.eq.11) open (5,file='STplot/plot11.dat')
               if (flplot.eq.12) open (5,file='STplot/plot12.dat')
               if (flplot.eq.13) open (5,file='STplot/plot13.dat')
               if (flplot.eq.14) open (5,file='STplot/plot14.dat')
               if (flplot.eq.15) open (5,file='STplot/plot15.dat')
               call minuit(fcn,chi2)
               close(5)
               if ((prob1.ge.chiref).and.(prob.le.chiref)) write(1,100) 
     .              Smin + 2*Sinc*(i - 1/2.d0 
     .                   + (prob - chiref)/(prob1 - prob)), Tpar
               if ((prob1.le.chiref).and.(prob.ge.chiref)) write(2,100) 
     .              Smin + 2*Sinc*(i - 1/2.d0 
     .                   + (prob - chiref)/(prob1 - prob)), Tpar
               prob1 = prob
 50         continue
 60      continue
      endif

 100  format(f8.4,f8.4)
      
      stop
      end
