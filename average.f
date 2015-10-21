      program sumrule

      implicit none
      real values(10), error(10),average,totaler
      integer i,imax
      logical corr

      imax = 10

c      data values/1.979,2.188,2.195,2.209,2.202,2.125,1.997,2.050,2.042,
c     .            2.075,1.946/

c      data error  /0.17, 0.15, 0.12, 0.14, 0.13, 0.11, 0.12, 0.11, 0.07,
c     .             0.08, 0.13/

c      data error  /0.18, 0.16, 0.13, 0.15, 0.14, 0.12, 0.13, 0.13, 0.09,
c     .             0.08, 0.14/

c      data values/291.029,290.6/
c      data error/0.780,1.1/

c      data values/20444.,20500./
c      data error/487.,227./

c      data values/21300.,22000.,20100.,20444./
c      data error/2700.,4800.,1900.,508./

c      data values/21300.,22000.,20100.,20490./
c      data error/2700.,4800.,1900.,251./

c      data values/-0.137, 0.052,-0.073,-0.064,-0.171,
c     .            -0.063,-0.004, 0.006, 0.112,-0.429/
c      data error/  0.250, 0.253, 0.146,0.0844, 0.266,
c     .             0.149, 0.146, 0.239, 0.212, 0.154/

      data values/-0.0045,-0.0233,-0.0046,-0.0116,-0.0125,
     .            -0.0114, 0.0019,-0.0015,-0.0051, 0.0040/
      data error/  0.0130, 0.0125, 0.0092, 0.0057, 0.0148,
     .             0.0095, 0.0097, 0.0142, 0.0124, 0.0105/

c      data values/ 291.183, 290.876/
c      data error/    0.979,   1.005/

      average = 0.
      totaler = 0.
      do 1 i = 1, imax
         average = average + values(i)/error(i)**2
         totaler = totaler + 1./error(i)**2
 1    continue
      average = average/totaler
      totaler = sqrt(1./totaler)
      
      print*,'average =',average,'+-',totaler

c      print*,(0.24*15. + 0.54*7.  + 1.06*4.  + 1.57*3.  + 1.32*2.  +
c     .        0.89*3.  + 1.03*4.  + 1.05+7.  + 0.39*15. + 0.15*30. +
c     .        0.24*40.)/137.**2/9/3.1416*2.

      stop
      end
