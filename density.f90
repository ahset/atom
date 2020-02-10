      subroutine density(nele,eigvec,nbasis,dens)
    
      implicit none
      integer nbasis, mi, ni, l, nele, i, j
      real*8 dens(nbasis,nbasis), eigvec(nbasis,nbasis)

      dens=0.0

! If many-electron atom and closed-shell

      if (nele.GT.1.and.MOD(nele,2)==0) then
        do mi=1,nbasis
          do ni=1,nbasis
            do l=1,(nele/2)
              dens(mi,ni)=dens(mi,ni)+2.0*(eigvec(mi,l)*eigvec(ni,l))
            enddo
          enddo
        enddo
      endif

! If many-electron atom and open-shell

      if (nele.GT.1.and.MOD(nele,2)==1) then
        do mi=1,nbasis
          do ni=1,nbasis
            do l=1,((nele-1)/2)
              dens(mi,ni)=dens(mi,ni)+2.0*(eigvec(mi,l)*eigvec(ni,l))+&
&eigvec(mi,1)*eigvec(ni,1)
            enddo
          enddo
        enddo
      endif

! If one-electron atom
    
      if (nele==1) then
        do mi=1,nbasis
          do ni=1,nbasis
            dens(mi,ni)=dens(mi,ni)+eigvec(mi,1)*eigvec(ni,1)
          enddo
        enddo
      endif

!dens(1,1) =0.008868
!dens(1,2) =0.030989
!dens(1,3) =0.064802
!dens(2,3) =0.226451
!dens(2,2) =0.108289
!dens(3,3) =0.473546
!do i=1,3
!    do j=1,3
!        dens(j,i)=dens(i,j)
!    enddo
!enddo

      return

      end subroutine density
