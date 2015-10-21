      subroutine s2plot
      
      implicit none
      double precision qq,shat,dgamma,alfaqq,msusy
      integer i
      include 'common.f'

      open (8,file='s2wplot/s2w_SM_2007.out',status='unknown')

      msusy = 250.d20

      do 10 i = 0,1000
         qq = 1.d4*dexp(-i/50.d0)
         call alfahat(qq,dgamma,alfaqq)
         write(8,100) qq,shat(qq)
 10   continue
      
 100  format(f14.7,f14.7)
      
      return
      end
