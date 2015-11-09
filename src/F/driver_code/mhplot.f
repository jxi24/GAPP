      program higgsplot

      implicit none
      integer i,j,hsteps
      double precision chi2,norm,probmh(0:1255),mhmin,mhinc,avprob
      external fcn
      external chi2
      include '../core/common.f'

      call mintio(7,8,9)

      open (5,file='higgs/mhplot_ICHEP08+FNAL.out',status='unknown')
      open (6,file='/dev/null',status='unknown')

      fwrite = .false.
      flprob = .true.
      flagmh = .false.
      flagmt = .true.
      flagmc = .true.
      flagal = .true.
      flagS  = .true.
      flagT  = .true.
      flgrho = .true.
      fkappa = .true.
      fsinth = .false.
      fzprim = .false.

      hsteps = 755
       mhmin =  99.4d0
       mhinc =   0.2d0
        norm =   0.0d0

      do 10 i = 1, hsteps 
         mh = mhmin + i*mhinc
         open (7,file='higgs/mh.dat',status='old')
         call minuit(fcn,chi2)
         close (7)
         probmh(i) = prob/mh
C        if (i.eq.1182) probmh(i) = 2.955d-17
         norm = norm + probmh(i)
 10   continue
      norm = norm/0.9986d0
      
      do 20 i = 1, hsteps
         if (i - 5*(i/5).eq.1) avprob = 0.d0
         avprob = avprob + probmh(i)
         if (i - 5*(i/5).eq.0) 
     .        write(5,*) mhmin + (i - 2)*mhinc, avprob/norm
 20   continue

      stop
      end
      
