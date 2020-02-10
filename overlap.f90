      subroutine overlap(nbasis, alpha, lap)
      
      implicit none
      integer nbasis, i, j
      real(8) alpha(nbasis), lap(nbasis,nbasis)
      
      do i=1,nbasis
        do j=1,nbasis
          lap(i,j)=(4.0*alpha(i)*alpha(j))**0.75/&
&(alpha(i)+alpha(j))**1.5
       end do
      end do
      
      return
      end subroutine overlap
