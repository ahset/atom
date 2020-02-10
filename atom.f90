      program atom
     
!     DFT calculation of 1 atom
!     Developed for Dr. Hélio Duarte's Density Functional
!    Theory class @ Chemistry Department - UFMG
!     By Taís Christofani

      implicit none
      integer allocation, nbasis, i, j,&
& zatom, nit,maxit, nele,charge
      real(8) ion,aff,chempot,trace, beta, etot, etot1, error, converg, exc, cxv,cxe,intvxcrho, t1, t2
      real(8), allocatable, DIMENSION(:,:) :: kin, lap, nuc,&
&eigvec, aux, dens, hcore, diag, dlap, xis, coul, fock,&
&gmini, vxc, dcore, dfock, hyb
      real(8), allocatable, DIMENSION(:) :: alpha, eigval, orb
      character(len=100)  :: title, output, doutput, poutput, fileplace, cwd

! Variables in the output:
! ZATOM > Atomic number
! NBASIS > Number of basis functions
! CHARGE > Charge of the atom
! ALPHA > Coefficients of the basis functions

! Open work files
      open (unit=1, file='input') ! Input information

! Read input information
      read(1,*) title
      read(1,*) zatom
      read(1,*) nele
      read(1,*) nbasis
      read(1,*) converg
      read(1,*) maxit
      read(1,*) cxv
      read(1,*) cxe
      read(1,*) beta   
      

      call getcwd(cwd)
      output = TRIM(title)//'.out'
      doutput = TRIM(title)//'.log'      
      poutput = TRIM(title)//'.plot'
      open (unit=2, file=TRIM(cwd)//'/results/'&
&//poutput) ! Radial density distribution data
      open (unit=3, file=TRIM(cwd)//'/results/'&
&//output) ! Output information 
      open (unit=4, file=TRIM(cwd)//'/results/'&
&//doutput) ! Detailed output information

      allocate(kin(nbasis,nbasis),alpha(nbasis),&
&lap(nbasis,nbasis),vxc(nbasis,nbasis), &
&nuc(nbasis,nbasis), eigval(nbasis), eigvec(nbasis,nbasis),&
&aux(nbasis,nbasis), hcore(nbasis,nbasis),dens(nbasis,nbasis),&
&diag(nbasis,nbasis), fock(nbasis,nbasis),&
&gmini(nbasis,nbasis),hyb(nbasis,nbasis),orb(nbasis),&
&dlap(nbasis,nbasis), xis(nbasis,nbasis), coul(nbasis,nbasis), &
&dcore(nbasis,nbasis),dfock(nbasis,nbasis), STAT=ALLOCATION)

      do i=1,nbasis
        read(1,*) alpha(i)
      end do


      write(3,*) "_______________________________________________________________"
      write(3,*) "|                                                             |"
      write(3,*) "| Program ATOM for calculation of the energy of a single atom |"
      write(3,*) "|     Developed as a part of  Dr. Hélio Duarte's DFT class    |"
      write(3,*) "|        ministered at the Chemistry Department of UFMG       |"
      write(3,*) "|                                                             |"
      write(3,*) "| By Taís Christofani                                         |"
      write(3,*) "|_____________________________________________________________|"
      write(3,*)  
      write(3,*)
      write(3,*) "******************"      
      write(3,*) "INPUT INFORMATION"
      write(3,*) "******************"
      write(3,*) 
      write(3,*) "Atomic number:"
      write(3,*) zatom
      write(3,*) 
      write(3,*) "Number of electron:"
      write(3,*) nele
      write(3,*) 
      write(3,*) "Number of basis functions:"
      write(3,*) nbasis
      write(3,*) 
      write(3,*) "Convergence criteria:"
      write(3,*) converg
      write(3,*) 
      write(3,*) "Maximum number of iterations:"
      write(3,*) maxit
      write(3,*)
      write(3,*) "Cx"
      write(3,*) cxe
      write(3,*)
      write(3,*) "Mixing coefficient for hybrid functional"
      write(3,*) beta
      write(3,*)
      write(3,*) "Gaussian functions coefficients:"
      do i=1,nbasis
        write(3,*) alpha(i)
      enddo
      write(3,*) 
      write(3,*)
      write(3,*) "******************"
      write(3,*) "OUTPUT INFORMATION"
      write(3,*) "******************"
      write(3,*)
      write(3,*) "* For detailed output information check .log file"
      write(3,*)

! Overlap integral

      call overlap(nbasis, alpha, lap)

      write(4,*)
      write(4,*) "Overlap matrix S:"
      do i=1,nbasis
        write(4,"(100g15.8)") (lap(i,j), j=1,nbasis)
      end do

! Nuclear integral

      call nuclear(nbasis, alpha, nuc, zatom)

      write(4,*)
      write(4,*) "Nuclear matrix N:"
      do i=1,nbasis
        write(4,"(100g15.8)") (nuc(i,j), j=1,nbasis)
      end do

! Kinect integral

      call kinect(nbasis, alpha, kin)

      write(4,*)
      write(4,*) "Kinect energy matrix:"
      do i=1,nbasis
        write(4,"(100g15.8)") (kin(i,j), j=1,nbasis)
      end do

! Obtain the transformation matrix X

      diag=lap
      call jacobi(diag,nbasis,eigval,eigvec)

      aux=0.0
      do i=1,nbasis
        aux(i,i)=eigval(i)**(-0.5)
      enddo

      xis=matmul(matmul(eigvec,aux),transpose(eigvec))

      write(4,*)
      write(4,*) "Transformation matrix X=S^(-1/2):"
      do i=1,nbasis
        write(4,"(100g15.8)") (xis(i,j), j=1,nbasis)
      end do

! Check if X is and orthogonalizing transformation matrix

      aux=matmul(matmul(xis,lap),xis)

      write(4,*)
      write(4,*) "XtSX:"
      do i=1,nbasis
        write(4,"(100g15.8)") (aux(i,j), j=1,nbasis)
      end do

! Obtaining a guess at the density matrix P

      hcore=kin+nuc

      write(4,*)
      write(4,*) "Initial Fock matrix F = Hcore:"
      do i=1,nbasis
        write(4,"(100g15.8)") (hcore(i,j), j=1,nbasis)
      end do

      call transfock(nbasis,hcore,eigvec,xis,eigval,orb)

! >>>>>>>>>>>>>>>>>>>>> SCF cycle <<<<<<<<<<<<<<<<<<<

      nit=1
      etot1=100.0
      error=100 

      write(4,*) "******************"
      write(4,*) "Entering SCF cycle"
      write(4,*) "******************"

      call cpu_time(t1)
      do while(nit<=maxit.and.error>converg)

        write(4,*)
        write(4,*) ":::::::::::::::: Step number",nit,"::::::::::::::::"

        call density(nele,eigvec,nbasis,dens)
        write(4,*)
        write(4,*) "Density matrix guess:"
        do i=1,nbasis
          write(4,"(100g15.8)") (dens(i,j), j=1,nbasis)
        end do
! Calculate the coulomb matrix
        
        call coulomb(nbasis,dens,alpha,coul)

        call hybrid(nbasis,dens,alpha,hyb)
        
        write(4,*)
        write(4,*) "Coulomb  matrix:"
        do i=1,nbasis
          write(4,"(100g15.8)") (coul(i,j), j=1,nbasis)
        end do

        write(4,*)
        write(4,*) "Exact exchange  matrix:"
        do i=1,nbasis
          write(4,"(100g15.8)") (hyb(i,j), j=1,nbasis)
        end do

! Calculate the exchange-correlation term

        call exco(dens,alpha,nbasis,cxv,cxe,vxc,exc)

        write(4,*)
        write(4,*) "Exchange-correlation potential matrix:"
        do i=1,nbasis
          write(4,"(100g15.8)") (vxc(i,j), j=1,nbasis)
        end do
! Calculate the fock matrix for closed-shell atoms

        if (nele.GT.1.and.MOD(nele,2)==0) then
          fock=hcore+coul+vxc+beta*(-0.5*hyb-vxc)
        endif

! Calculate the fock matrix for the hydrogen atom

        if (nele==1) then
          fock=hcore+coul+vxc+beta*(-hyb-vxc)
        endif

        write(4,*)
        write(4,*) "Fock matrix F:"
        do i=1,nbasis
          write(4,"(100g15.8)") (fock(i,j), j=1,nbasis)
        end do

        call transfock(nbasis,fock,eigvec,xis,eigval,orb)

! Calculate the energy at this step
        
        call energy(nbasis,hcore,kin,coul,hyb,dens,exc,etot,beta,nele)
        nit=nit+1
          error=abs(etot1-etot)
          etot1=etot
      enddo
      call cpu_time(t2)
 
! >>>>>>>>>>>>>>>>>>>  End of SCF <<<<<<<<<<<<<<<<


      write(4,*)
      write(4,*) "******************"
      write(4,*) "End of  SCF cycle"
      write(4,*) "******************"

      write(3,*) "Number of SCF cycles for convergence:"
      write(3,*) 
      write(3,*) nit-1


! Mulliken population analysis

      charge=0.0
      trace=0.0

      aux=matmul(dens,lap)
      do i=1,nbasis
        trace=trace+aux(i,i)
      enddo

      charge=(zatom-trace)

      write(3,*)
      write(3,*) "Charge:"
      write(3,*) charge

      if (nele.GT.1.and.MOD(nele,2)==0) then
        ion=-orb(nele/2)
        aff=-orb(nele/2+1)
      endif

      if (nele==1) then
        ion=-orb(nele)
        aff=-orb(nele+1)
      endif

      chempot=-(ion+aff)/2.0

      write(3,*)
      write(3,*) "Ionization potential:"
      write(3,*) ion
      write(3,*)
      write(3,*) "Electron affinity:"
      write(3,*) aff
      write(3,*)
      write(3,*) "Chemical potential"
      write(3,*) chempot
      write(3,*) 
      write(3,*) "Density:"
      do i=1,nbasis
        write(3,"(100g15.8)") (dens(i,j), j=1,nbasis)
      enddo
      write(3,*)
      write(3,*) "Orbital energies:"
      do i=1,nbasis
        write(3,"(100g15.8)") i,"",orb(i)
      enddo
      write(3,*)
      write(3,*) "Total energy:"
      write(3,*) etot


      write(3,*)
      write(3,*) "---------------------------------------------"
      write(3,*) "CPU time:",t2-t1,"seconds"
      
      close(1)
      close(3)
      close(4)
     
      deallocate (dlap, gmini, dens, alpha, kin,&
& lap, nuc, fock, eigval, eigvec, aux, hcore, diag, xis,&
& coul, orb,vxc, dcore, dfock, hyb, STAT=ALLOCATION)

      end program
