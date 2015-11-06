      program smfit

      implicit none
      double precision chi2
      external fcn
      external chi2
      include 'common.f'

CC--K      call mintio(5,7,9)

      open (5,file='smfit.dat',status='old')
      open (6,file='/dev/null',status='unknown')
      open (7,file='uu-d.out',status='unknown')
      open (50,file='uu-d.pull',status='unknown')
      open (51,file='summary.out',access = 'append',status='unknown')

C Kai023 ---------------------------------
C  Suppress unimportant output and allow to read in NP parameters from
C  the input file.
        
C     flgfitx = .false.
C     flgtph  = .false.
C     flgs2b  = .false.

C     fitx =   400.0d0
C     fittph = 4.26d0
C     fits2b = 0.000d0

      flgfitx = .true.
      flgtph  = .true.
      flgs2b  = .true.

      flfout = .true.

C     fwrite = .true.
      fwrite = .false.
C ----------------------------------------

      fsplot = .true.
      flprob = .true.
C ----------------------------------------
      flagmh = .true.
C     flagmh = .false.
C     mh = 25.0d0
C ----------------------------------------
      flagmt = .true.
      flagmc = .true.
      flagal = .true.
      flagS  = .true.
      flagT  = .true.
      flgrho = .true.
      fkappa = .true.
      fsinth = .true.
      fzprim = .true.

      sigma = 5.d0

CC--K      call minuit(fcn,chi2)

      stop
      end
