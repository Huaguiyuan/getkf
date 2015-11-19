!
!  This simple program calculates accurate Fermi vectors from Wannier output
!    A moderate fs.bxsf is required (produced by wannier90)
!     also needs the wannier90_hr.dat file (Hamiltonian produced by wannier90)
!    Then using the mesh in fs.bxsf, accurate k_F^y is calculated.
!
!  Usage: getkf.x fs.bxsf wannier90_hr.dat offset
!
PROGRAM getkf
  !
  USE constants, ONLY : dp, eps5
  USE banddata,  ONLY : read_band_from_bxsf, ef, eig, kvec, &
                        nbnd, nkpt, nkx, nky, nkz
  USE wanndata,  ONLY : read_ham, finalize_wann
  !
  IMPLICIT NONE
  !
  character(len=80) arg
  !
  integer nn ! Band offset
  integer ix, iy, iz, ib, ii, jj
  real(dp) kf(1:3)
  !
  CALL getarg(1, arg)
  !
  CALL read_band_from_bxsf(arg)
  !
  CALL getarg(2, arg)
  !
  CALL read_ham(arg)
  !
  CALL getarg(3, arg)
  !
  read(arg, *) nn
  !
  do ib=1, nbnd
    do ii=1, nkpt
      eig(ii, ib)=eig(ii, ib)-ef
    enddo
  enddo
  !
  do ib=1, nbnd
    do iz=0, nkz-1
      do ix=0, nkx-1
        do iy=0, nky-2
          ii=ix*nky*nkz+iy*nkz+iz+1
          jj=ix*nky*nkz+(iy+1)*nkz+iz+1
          if ( eig(ii, ib)*eig(jj, ib) <= 0 ) then
            CALL solve_kf(kf, kvec(ii, :), kvec(jj, :), ib+nn)
!            kf(:)=(eig(ii,ib)*kvec(jj,:)-eig(jj, ib)*kvec(ii, :))/(eig(ii, ib)-eig(jj, ib))
!            CALL solve_kf(kf, ef, ib+nn, eps5)
            write(*,'(1I4,3F16.9)') ib, kf(:)
          endif
        enddo ! iy
      enddo ! ix
    enddo ! iz
  enddo ! ib
  !
  CALL finalize_wann
  !
END PROGRAM
