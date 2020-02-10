      subroutine energy(nbasis,hcore,kin,coul,hyb,dens,exc,etot,beta,nele)
    
      implicit none
      integer nbasis, i, j, nele
      real(8) coul(nbasis,nbasis), hcore(nbasis,nbasis), kin(nbasis,nbasis),&
& dens(nbasis,nbasis), hyb(nbasis,nbasis), ecore, exchf, ecoul, exc, etot, enuc, ekin, beta

      ecore=0.0
      ecoul=0.0
      enuc=0.0
      ekin=0.0

      if (nele.GT.1.and.MOD(nele,2)==0) then
        do i=1,nbasis
          do j=1,nbasis
            ekin=ekin+dens(j,i)*kin(i,j)
            ecore=ecore+dens(j,i)*hcore(i,j)
            ecoul=ecoul+dens(j,i)*0.5*coul(i,j)
            exchf=exchf+dens(j,i)*0.25*hyb(i,j)
          enddo
        enddo    
        enuc=ecore-ekin
        etot=ecore+ecoul+exc+beta*(-exchf-exc)
      endif

      if (nele==1) then
        do i=1,nbasis
          do j=1,nbasis
            ekin=ekin+dens(j,i)*kin(i,j)
            ecore=ecore+dens(j,i)*hcore(i,j)
            ecoul=ecoul+dens(j,i)*0.5*coul(i,j)
            exchf=exchf+dens(j,i)*0.5*hyb(i,j)
          enddo
        enddo
        enuc=ecore-ekin
        etot=ecore+ecoul+exc+beta*(-exchf-exc)
      endif


      write(4,*)
      write(4,*) "Kinetic energy:"
      write(4,*) ekin
      write(4,*)
      write(4,*) "Nuclear-elctron potential energy:"
      write(4,*) enuc
      write(4,*)
      write(4,*) "===> One-electron energy:"
      write(4,*) ecore
      write(4,*)
      write(4,*) "Coulomb energy:"
      write(4,*) ecoul
      write(4,*)
      write(4,*) "Exchange-correlation energy"
      write(4,*) exc
      write(4,*)
      write(4,*) "HF exact exchange energy:"
      write(4,*) exchf
      write(4,*)
      write(4,*) "===> Two-electron energy:"
      write(4,*) exc+ecoul
      write(4,*)
      write(4,*) "===> Total energy:"
      write(4,*) etot
      return
      end subroutine energy
