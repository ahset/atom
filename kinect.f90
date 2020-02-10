      subroutine kinect(nbasis, alpha, kin)
    
      implicit none
      integer nbasis, i, j
      real*8 alpha(nbasis), kin(nbasis,nbasis)
      do i=1,nbasis
        do j=1,nbasis
          kin(i,j)=6.0*2.0**0.5*alpha(i)**0.75*alpha(j)**&
&1.75*(1.0-alpha(j)/(alpha(i)+alpha(j)))/(alpha(i)+alpha(j))**1.5
        end do
      end do

      return

      end subroutine kinect
