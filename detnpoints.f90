      subroutine detnpoints(stepsize,dens,alpha,nbasis,npoints,rrmat,rhomat,rhorad)

      integer npoints, i, j, k, l, nbasis
      real(8) intprec, stepsize, phia, phib, alpha(nbasis)
      real(8) rho,rr,pi,dens(nbasis,nbasis)
      real(8)  rrmat(1000),rhomat(1000), rhorad(1000)

      pi=4.0*atan(1.0)
      npoints=0

      rr=0.0
      rho=10.0  
      intprec=10E-20
      stepsize=0.05

      do while (rho.GT.intprec.and.npoints.LT.500)
      rho=0.0
      npoints=npoints+1
      l=npoints
        do i=1,nbasis
          phia=(2.0*alpha(i)/pi)**0.75*exp(-alpha(i)*rr**2.0)
          do j=1,nbasis
            phib=(2.0*alpha(j)/pi)**0.75*exp(-alpha(j)*rr**2.0)
            rho=rho+dens(i,j)*phia*phib
          enddo
        enddo
        rrmat(l)=rr
        rhomat(l)=rho
        rhorad(l)=rho*4.0*pi*rr**2.0
        rr=rr+stepsize
      enddo
   
      write (2,*) "r vs rho"
      write (2,*)
      do i=1,npoints
        write (2,*) rrmat(i),rhomat(i)*4.0*pi*rrmat(i)**2.0
      enddo
 
      return

      end subroutine detnpoints
