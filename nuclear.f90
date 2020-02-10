      subroutine nuclear(nbasis, alpha, nuc, zatom)
      implicit none
      integer nbasis, i, j, zatom
      real(8) alpha(nbasis), nuc(nbasis,nbasis)
      real pi
      pi = 4.0*atan(1.0)
      do i=1,nbasis
        do j=1,nbasis
          nuc(i,j)=-2.0**2.5*pi**(-0.5)*zatom*(alpha(i)*&
&alpha(j))**0.75/(alpha(i)+alpha(j))
       end do
      end do
      
      return
      end subroutine nuclear
