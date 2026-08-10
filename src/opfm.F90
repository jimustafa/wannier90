!-*- mode: F90 -*-!
!------------------------------------------------------------!
! Copyright (C) 2026 Wannier Developer Group                 !
!                                                            !
! This library is free software; you can redistribute it     !
! and/or modify it under the terms of the GNU Lesser General !
! Public License as published by the Free Software           !
! Foundation; either version 2.1 of the License, or (at your !
! option) any later version.                                 !
!                                                            !
! This library is distributed in the hope that it will be    !
! useful,but WITHOUT ANY WARRANTY; without even the implied  !
! warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR    !
! PURPOSE.  See the GNU Lesser General Public License for    !
! more details.                                              !
!                                                            !
! You should have received a copy of the GNU Lesser General  !
! Public License along with this library; if not, see        !
! <https://www.gnu.org/licenses/>.                           !
!                                                            !
! The webpage of the Wannier90 code is                       !
! <https://www.wannier.org>.                                 !
!                                                            !
! The Wannier90 code is hosted on GitHub                     !
! <https://github.com/wannier-developers/wannier90>          !
!------------------------------------------------------------!
!                                                            !
!  w90_opfm: Optimized Projection Functions Method           !
!                                                            !
!  Refs: 10.1103/PhysRevB.92.165134                          !
!        10.1103/PhysRevB.94.125151                          !
!                                                            !
!------------------------------------------------------------!

module w90_opfm

  use w90_constants, only: dp, cmplx_0, cmplx_1, cmplx_i, pi
  use w90_error, only: w90_error_type, set_error_alloc, set_error_dealloc, set_error_fatal
  use w90_comms, only: w90_comm_type, comms_allreduce, comms_bcast, mpirank
  use w90_types, only: kmesh_info_type

  implicit none

  private

  public :: opfm_lowdin
  public :: opfm_setup
  public :: opfm_random_unitary

contains

  subroutine opfm_random_unitary(n, w, error, comm)
    !================================================!
    !
    !! Generate a Haar-random n x n unitary matrix via QR decomposition of a
    !! complex Ginibre matrix, with the R-diagonal phase correction needed
    !! for the Q factor to be Haar-distributed. Generated on the root rank
    !! and broadcast so all ranks agree.
    !
    !================================================!

    integer, intent(in) :: n
    complex(kind=dp), intent(out) :: w(:, :)
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    complex(kind=dp), allocatable :: z(:, :), cwork(:), tau(:)
    real(kind=dp), allocatable :: u1(:, :), u2(:, :)
    complex(kind=dp) :: cwork_query(1)
    integer :: i, info, ierr, lwork

    if (mpirank(comm) == 0) then
      allocate (z(n, n), u1(n, n), u2(n, n), tau(n), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating workspace in opfm_random_unitary', comm)
        return
      end if

      call random_number(u1)
      call random_number(u2)
      z = sqrt(-2.0_dp*log(1.0_dp - u1))*cos(2.0_dp*pi*u2)/sqrt(2.0_dp) + &
          cmplx_i*sqrt(-2.0_dp*log(1.0_dp - u1))*sin(2.0_dp*pi*u2)/sqrt(2.0_dp)

      call zgeqrf(n, n, z, n, tau, cwork_query, -1, info)
      lwork = max(1, int(real(cwork_query(1))))
      allocate (cwork(lwork), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error allocating cwork in opfm_random_unitary', comm)
        return
      end if
      call zgeqrf(n, n, z, n, tau, cwork, lwork, info)
      if (info /= 0) then
        call set_error_fatal(error, 'Error in ZGEQRF in opfm_random_unitary', comm)
        return
      end if

      do i = 1, n
        w(i, i) = z(i, i)
      end do

      call zungqr(n, n, n, z, n, tau, cwork, lwork, info)
      if (info /= 0) then
        call set_error_fatal(error, 'Error in ZUNGQR in opfm_random_unitary', comm)
        return
      end if

      do i = 1, n
        w(:, i) = z(:, i)*w(i, i)/abs(w(i, i))
      end do

      deallocate (z, u1, u2, tau, cwork, stat=ierr)
      if (ierr /= 0) then
        call set_error_dealloc(error, 'Error deallocating workspace in opfm_random_unitary', comm)
        return
      end if
    end if

    call comms_bcast(w(1, 1), n*n, error, comm)
    if (allocated(error)) return
  end subroutine opfm_random_unitary

  subroutine opfm_lowdin(a, u, error, comm)
    !================================================!
    !
    !! Lowdin-orthonormalize each a(:,:,k) (n_rows x n_cols, n_rows <= n_cols):
    !! u(k) = Z(k).V(k)^dagger from the economy SVD a(k) = Z(k).Sigma(k).V(k)^dagger
    !
    !================================================!

    complex(kind=dp), intent(in)  :: a(:, :, :)
    complex(kind=dp), intent(out) :: u(:, :, :)
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer :: n_rows, n_cols, num_kpts, ik, info, ierr, lwork
    real(kind=dp), allocatable :: svals(:), rwork(:)
    complex(kind=dp), allocatable :: cwork(:), cz(:, :), cvdag(:, :)
    complex(kind=dp) :: cwork_query(1)

    n_rows = size(a, 1)
    n_cols = size(a, 2)
    num_kpts = size(a, 3)

    allocate (svals(n_rows), rwork(5*n_rows), cz(n_rows, n_rows), cvdag(n_rows, n_cols), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating workspace in opfm_lowdin', comm)
      return
    end if

    call zgesvd('S', 'S', n_rows, n_cols, a(1, 1, 1), n_rows, svals, cz, n_rows, cvdag, n_rows, &
                cwork_query, -1, rwork, info)
    lwork = max(1, int(real(cwork_query(1))))
    allocate (cwork(lwork), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating cwork in opfm_lowdin', comm)
      return
    end if

    do ik = 1, num_kpts
      u(:, :, ik) = a(:, :, ik)
      call zgesvd('S', 'S', n_rows, n_cols, u(1, 1, ik), n_rows, svals, cz, n_rows, cvdag, n_rows, &
                  cwork, lwork, rwork, info)
      if (info /= 0) then
        call set_error_fatal(error, 'Error in ZGESVD in opfm_lowdin', comm)
        return
      end if
      call zgemm('N', 'N', n_rows, n_cols, n_rows, cmplx_1, cz, n_rows, cvdag, n_rows, cmplx_0, &
                 u(:, :, ik), n_rows)
    end do
  end subroutine opfm_lowdin

  subroutine opfm_setup(kmesh_info, a_matrix, m_matrix_loc, lambda, include_bweights, nkrank, &
                        global_k, num_wann, num_proj, num_kpts, sx_matrix, mx_matrix, ww, &
                        error, comm)
    !================================================!
    !
    !! Build the codiagonalization inputs for the general (num_proj >= num_wann)
    !! case:
    !!   u_a(:,:,k)         = Lowdin(a_matrix(:,:,k))
    !!   sx_matrix(:,:,k)   = sqrt(lambda) . (a_matrix(:,:,k)^dagger . a_matrix(:,:,k) - I)  (num_proj x num_proj)
    !!   mx_matrix(:,:,idx) = sqrt(wb) . u_a(:,:,k)^dagger . m_matrix_loc . u_a(:,:,k+b)
    !!   ww                 = identity (num_proj x num_proj)
    !! mx_matrix and sx_matrix are kept separate (not concatenated) since the
    !! codiagonalization sweep must maximize diag2(mx_matrix) while minimizing
    !! diag2(sx_matrix).
    !
    !================================================!

    type(kmesh_info_type), intent(in) :: kmesh_info
    complex(kind=dp), intent(in)  :: a_matrix(:, :, :)
    complex(kind=dp), intent(in)  :: m_matrix_loc(:, :, :, :)
    real(kind=dp), intent(in)     :: lambda
    logical, intent(in)           :: include_bweights
    integer, intent(in)           :: nkrank
    integer, intent(in)           :: global_k(:)
    integer, intent(in)           :: num_wann
    integer, intent(in)           :: num_proj
    integer, intent(in)           :: num_kpts
    complex(kind=dp), intent(out) :: sx_matrix(:, :, :)
    complex(kind=dp), intent(out) :: mx_matrix(:, :, :)
    complex(kind=dp), intent(out) :: ww(:, :)
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    complex(kind=dp), allocatable :: cvdag(:, :), u_a(:, :, :)
    real(kind=dp) :: bweight, lambda_eff
    integer :: nkp, nkp2, nkp_loc, nn, idx, ik, i, ierr

    allocate (cvdag(num_proj, num_wann), u_a(num_wann, num_proj, num_kpts), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error allocating workspace in opfm_setup', comm)
      return
    end if

    call opfm_lowdin(a_matrix, u_a, error, comm)
    if (allocated(error)) return

    ! when the m-stack slices are b-weighted, rescale lambda by sum(wb) so the
    ! S-penalty stays on the same footing as the (now b-weighted) m-term
    lambda_eff = lambda
    if (include_bweights) lambda_eff = lambda*sum(kmesh_info%wb(1:kmesh_info%nntot))

    do ik = 1, num_kpts
      call zgemm('C', 'N', num_proj, num_proj, num_wann, cmplx_1, a_matrix(:, :, ik), num_wann, &
                 a_matrix(:, :, ik), num_wann, cmplx_0, sx_matrix(:, :, ik), num_proj)
      do i = 1, num_proj
        sx_matrix(i, i, ik) = sx_matrix(i, i, ik) - cmplx_1
      end do
      sx_matrix(:, :, ik) = sqrt(lambda_eff)*sx_matrix(:, :, ik)
    end do

    mx_matrix = cmplx_0
    do nkp_loc = 1, nkrank
      nkp = global_k(nkp_loc)
      do nn = 1, kmesh_info%nntot
        nkp2 = kmesh_info%nnlist(nkp, nn)
        idx = (nkp - 1)*kmesh_info%nntot + nn
        bweight = 1.0_dp
        if (include_bweights) bweight = sqrt(kmesh_info%wb(nn))
        call zgemm('C', 'N', num_proj, num_wann, num_wann, cmplx_1, u_a(:, :, nkp), num_wann, &
                   m_matrix_loc(:, :, nn, nkp_loc), num_wann, cmplx_0, cvdag, num_proj)
        call zgemm('N', 'N', num_proj, num_proj, num_wann, bweight*cmplx_1, cvdag, &
                   num_proj, u_a(:, :, nkp2), num_wann, cmplx_0, mx_matrix(:, :, idx), num_proj)
      end do
    end do

    call comms_allreduce(mx_matrix(1, 1, 1), num_proj*num_proj*num_kpts*kmesh_info%nntot, 'SUM', error, comm)
    if (allocated(error)) return

    ww = cmplx_0
    do i = 1, num_proj
      ww(i, i) = cmplx_1
    end do

    deallocate (cvdag, u_a, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error deallocating workspace in opfm_setup', comm)
      return
    end if
  end subroutine opfm_setup

end module w90_opfm
