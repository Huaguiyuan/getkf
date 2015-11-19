!
!   wanndata.f90
!   
!
!   Created by Chao Cao on 01/03/14.
!   Copyright 2013 __MyCompanyName__. All rights reserved.
!

#include 'lapack.f90'

MODULE wanndata
  !
  use constants
  !
  IMPLICIT NONE

  INTEGER norb
  INTEGER nrpt

  COMPLEX(DP), ALLOCATABLE :: ham(:,:,:)
  REAL(DP), ALLOCATABLE :: weight(:)
  REAL(DP), ALLOCATABLE :: rvec(:,:)
  !
CONTAINS

SUBROUTINE read_ham(fn)
!
  USE constants
  !
  IMPLICIT NONE
  !
  CHARACTER(len=80) fn
  INTEGER irpt, iorb, jorb, t1, t2, t3, t4, t5
  INTEGER, ALLOCATABLE :: wt(:)
  REAL(DP) a, b
  !
  open(unit=fin, file=fn)
  !
  read(fin, *)
  read(fin, *) norb
  read(fin, *) nrpt
  !
  write(stdout, *) " #  Dimensions:"
  write(stdout, *) "    # of orbitals:", norb
  write(stdout, *) "    # of real-space grid:", nrpt
  !
  allocate(ham(1:norb, 1:norb, 1:nrpt))
  allocate(weight(1:nrpt))
  allocate(rvec(1:3, 1:nrpt))
  !
  allocate(wt(1:nrpt))
  read(fin, '(15I5)') (wt(irpt),irpt=1,nrpt)
  weight(:)=wt(:)
  deallocate(wt)
  !
  do irpt=1, nrpt
    do iorb=1, norb
      do jorb=1, norb
        read(fin, *) t1, t2, t3, t4, t5, a, b
        if ((iorb.eq.1).and.(jorb.eq.1)) then
          rvec(1, irpt)=t1
          rvec(2, irpt)=t2
          rvec(3, irpt)=t3
        endif
        ham(iorb, jorb, irpt)=CMPLX(a,b)
      enddo
    enddo
  enddo
  !
  close(unit=fin)
  !
END SUBROUTINE

SUBROUTINE finalize_wann()
  !
  IMPLICIT NONE
  !
  if (allocated(ham)) deallocate(ham)
  if (allocated(weight)) deallocate(weight)
  if (allocated(rvec)) deallocate(rvec)
  !
END SUBROUTINE

SUBROUTINE interpolate_band(ek, kpt)
  !
  USE constants, ONLY: dp, twopi, cmplx_i, cmplx_0
  USE lapack95,  ONLY: heev
  !
  IMPLICIT NONE
  !
  real(dp),dimension(1:norb):: ek
  real(dp),dimension(1:3):: kpt
  !
  real(dp) rdotk
  complex(dp) fact
  complex(dp), allocatable :: work(:, :)
  !
  integer ir, info
  !
  allocate(work(1:norb, 1:norb))
  work(:,:)=cmplx_0
  !
  do ir=1, nrpt
    rdotk=SUM(kpt(:)*rvec(:, ir))
    fact=exp(-cmplx_i*twopi*rdotk)/weight(ir)
    work(:,:)=work(:,:)+ham(:,:,ir)*fact
  enddo ! ir
  !
  call heev(work, ek, 'V', 'U', info)
  !
  deallocate(work)
  !
END SUBROUTINE

END MODULE

