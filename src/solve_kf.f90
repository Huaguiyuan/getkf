SUBROUTINE solve_kf(kf, kv1, kv2, ib)
  !
  USE constants, ONLY: dp, eps9
  USE wanndata,  ONLY: interpolate_band, norb
  USE banddata,  ONLY: ef
  !
  IMPLICIT NONE
  !
  real(dp),dimension(1:3) :: kf, kv1, kv2
  integer   ib
  !
  real(dp),dimension(1:3) :: k1, k2
  real(dp)  e1, e2, eig
  real(dp),allocatable :: ek(:)
  integer ii
  !
  allocate(ek(1:norb))
  !
  k1=kv1
  k2=kv2
  CALL interpolate_band(ek, k1)
  e1=ek(ib)-ef
  CALL interpolate_band(ek, k2)
  e2=ek(ib)-ef
  eig=1.d0
  !
  do while(abs(eig)>eps9)
    kf(:)=(e1*k2(:)-e2*k1(:))/(e1-e2)
    CALL interpolate_band(ek, kf)
    eig=ek(ib)-ef
    if (eig*e1<=0) then
      k2=kf
      e2=eig
    else
      k1=kf
      e1=eig
    endif
  enddo
  !
  deallocate(ek)
  !
END SUBROUTINE
