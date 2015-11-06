      program sort

      implicit none
      integer i,N(200)
      double precision X(200),Y(200)

      open (1,file='mncont.out',status='old')
      open (2,file='npplots/ST_all.out',status='unknown')

       read(1,*)   (N(i),X(i),Y(i),N(i+10),X(i+10),Y(i+10),i=1, 10)
      write(2,100) (     X(i),Y(i)                        ,i=1, 20)

 100  format(X,f10.4,f10.4)

      stop
      end
