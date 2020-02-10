      subroutine coulomb(nbasis, dens, alpha, coul)
    
      implicit none
      integer nbasis, sigma, lambda, mi, ni
      real(8) alpha(nbasis), coul(nbasis,nbasis)
      real(8) dens(nbasis,nbasis), pi, som

      pi = 4.0*atan(1.0)
   
      do mi=1,nbasis
        do ni=1,nbasis
          som=0.0
          do sigma=1,nbasis
            do lambda=1,nbasis
              som=som+dens(sigma,lambda)*16.0*pi**(-0.5)*&
&(alpha(mi)*alpha(ni)*alpha(sigma)*alpha(lambda))**(0.75)/&
&((alpha(sigma)+alpha(lambda))*(alpha(ni)+alpha(mi))*&
&(alpha(mi)+alpha(ni)+alpha(sigma)+alpha(lambda))**(0.5))
            enddo
          enddo
          coul(mi,ni)=som
        enddo
      enddo
      
      return

      end subroutine coulomb
