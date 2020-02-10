      subroutine transfock(nbasis,fock,eigvec,xis,eigval,orb)
    
      implicit none
      integer nbasis, i, j
      real(8) fock(nbasis,nbasis), xis(nbasis,nbasis),&
&aux(nbasis,nbasis), eigvec(nbasis,nbasis), diag(nbasis,nbasis),&
&eigval(nbasis),orb(nbasis)

! Calculate the transformed fock matrix F'

      aux=matmul(matmul(transpose(xis),fock),xis)

      write(4,*)
      write(4,*) "Transformed Fock matrix F':"
      do i=1,nbasis
        write(4,"(100g15.8)") (aux(i,j), j=1,nbasis)
      end do

! Diagonalize the transformed fock matrix

      diag=aux

      call jacobi(diag, nbasis, eigval, eigvec)

      orb=eigval
   
      write(4,*)
      write(4,*) "Orbitals energies:"
      do i=1,nbasis
        write(4,"(100g15.8)") i,eigval(i)
      end do
 
      write(4,*)
      write(4,*) "Transformed coefficients C':"
      do i=1,nbasis
        write(4,"(100g15.8)") (eigvec(i,j), j=1,nbasis)
      end do

! Calculate the new coefficients matrix

      eigvec=matmul(xis,eigvec)

      write(4,*)
      write(4,*) "New coefficients C=XC':"
      do i=1,nbasis
        write(4,"(100g15.8)") (eigvec(i,j), j=1,nbasis)
      end do

      return
      end subroutine transfock
