MODULE banddata
  !
  USE constants, ONLY : dp
  !
  IMPLICIT NONE
  !
  integer nkpt, nbnd
  integer nkx, nky, nkz
  real(dp) ef
  real(dp),allocatable :: kvec(:,:)
  real(dp),allocatable :: eig(:,:)
  !
 CONTAINS
 !
 SUBROUTINE create_mesh()
  !
  USE constants, ONLY : dp
  !
  IMPLICIT NONE
  !
  integer ix, iy, iz
  !
  nkpt=nkx*nky*nkz
  allocate(kvec(nkpt, 3))
  do ix=0, nkx-1
    do iy=0, nky-1
      do iz=0, nkz-1
        kvec(ix*nky*nkz+iy*nkz+iz+1,1)=(ix*1.d0)/(nkx-1)
        kvec(ix*nky*nkz+iy*nkz+iz+1,2)=(iy*1.d0)/(nky-1)
        kvec(ix*nky*nkz+iy*nkz+iz+1,3)=(iz*1.d0)/(nkz-1)
      enddo
    enddo
  enddo
  !
 END SUBROUTINE
 !
 SUBROUTINE read_band_from_bxsf(fn)
  !
  USE constants, ONLY : dp, fin
  !
  IMPLICIT NONE
  !
  character(len=80):: fn
  character(len=80):: dummy
  integer ii, ik, ib
  !
  open(unit=fin, file=fn)
  !
  do ii=1, 8
    read(fin, *) dummy
  enddo
  !
  read(fin, '(1A20, 1F)') dummy, ef
  !
  write(*,*) "# Ef=",ef
  !
  do ii=1, 5
    read(fin, *)
  enddo
  !
  read(fin, *) nbnd
  read(fin, *) nkx, nky, nkz
  CALL create_mesh()
  !
  do ii=1, 4
    read(fin, *)
  enddo
  !
  allocate(eig(nkpt, nbnd))
  !
  do ib=1, nbnd
    read(fin, *)
    do ik=1, nkpt
      read(fin, *) eig(ik, ib)
    enddo
  enddo
  !
  close(unit=fin)
  !
 END SUBROUTINE
 !
 SUBROUTINE destroy_banddata
  !
  IMPLICIT NONE
  !
  deallocate(eig, kvec)
  !
 END SUBROUTINE
 !
END MODULE
