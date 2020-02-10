      subroutine exco(dens,alpha,nbasis,cxv,cxe,vxc,exc)

      integer npoints, i, j, k, nbasis
      real(8) integral, stepsize, phia, phib, alpha(nbasis),cxv,cxe
      real(8) RR,rho,pi,dens(nbasis,nbasis), som, exc
      real(8) intxc(1000),vxc(nbasis,nbasis),rhomat(1000),rrmat(1000),rhorad(1000)

      pi=4.0*atan(1.0)

      call detnpoints(stepsize,dens,alpha,nbasis,npoints,rrmat,rhomat,rhorad)

! Calculate POTxc
      integral=0.0
      do i=1,nbasis
        do j=1,nbasis
          do k=1,npoints
            phia=(2.0*alpha(i)/pi)**0.75*exp(-alpha(i)*rrmat(k)**2.0)
            phib=(2.0*alpha(j)/pi)**0.75*exp(-alpha(j)*rrmat(k)**2.0)
            intxc(k)=-cxv*phia*rhomat(k)**(1.0/3.0)*phib*4.0*pi*rrmat(k)**2.0
          enddo
          call simpson(intxc,npoints,stepsize,integral)
          vxc(i,j)=integral
        enddo
      enddo

! Calculate Exc
      integral=0.0
      do k=1,npoints
        intxc(k)=-cxe*rhomat(k)**(4.0/3.0)*4.0*pi*rrmat(k)**2.0
      enddo
      call simpson(intxc,npoints,stepsize,integral)
      exc=integral

! Integrating RHO to give electron number    
      call simpson(rhorad,NPOINTS,STEPSIZE,integral)
      write(4,*) "Calculated number of electrons:"
      write(4,*) integral
            
      return
      end subroutine exco
