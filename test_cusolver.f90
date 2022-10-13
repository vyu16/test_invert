SUBROUTINE test_cusolver_real(ndim,mat)
   !
   USE cudafor
   USE cusolverdn
   !
   IMPLICIT NONE
   !
   INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14,200)
   !
   INTEGER, INTENT(IN) :: ndim
   REAL(DP), INTENT(INOUT) :: mat(ndim,ndim)
   !
   INTEGER :: lwork
   INTEGER :: ierr
   INTEGER :: cr
   INTEGER :: t1
   INTEGER :: t2
   INTEGER, ALLOCATABLE :: ipiv(:)
   INTEGER, DEVICE :: info_d
   INTEGER, DEVICE, ALLOCATABLE :: ipiv_d(:)
   !
   REAL(DP), ALLOCATABLE :: work(:)
   REAL(DP), DEVICE, ALLOCATABLE :: mat_d(:,:)
   REAL(DP), DEVICE, ALLOCATABLE :: work_d(:)
   !
   TYPE(cusolverDnHandle) :: cusolver_h
   !
   ierr = cusolverDnCreate(cusolver_h)
   ierr = cusolverDnDgetrf_bufferSize(cusolver_h,ndim,ndim,mat_d,ndim,lwork)
   !
   ALLOCATE(ipiv(ndim))
   ALLOCATE(work(1))
   ALLOCATE(mat_d(ndim,ndim))
   ALLOCATE(ipiv_d(ndim))
   ALLOCATE(work_d(lwork))
   !
   mat_d(:,:) = mat
   !
   ierr = cudaDeviceSynchronize()
   !
   CALL SYSTEM_CLOCK(COUNT_RATE=cr)
   CALL SYSTEM_CLOCK(t1)
   !
   ierr = cusolverDnDgetrf(cusolver_h,ndim,ndim,mat_d,ndim,work_d,ipiv_d,info_d)
   !
   CALL DGETRI(ndim,mat,ndim,ipiv,work,-1,ierr)
   !
   lwork = CEILING(work(1))
   !
   DEALLOCATE(work)
   ALLOCATE(work(lwork))
   !
   ipiv(:) = ipiv_d
   mat(:,:) = mat_d
   !
   CALL DGETRI(ndim,mat,ndim,ipiv,work,lwork,ierr)
   !
   CALL SYSTEM_CLOCK(t2)
   !
   PRINT *,"time2:",REAL(t2-t1)/REAL(cr)
   !
   DEALLOCATE(ipiv)
   DEALLOCATE(work)
   DEALLOCATE(mat_d)
   DEALLOCATE(ipiv_d)
   DEALLOCATE(work_d)
   !
   ierr = cusolverDnDestroy(cusolver_h)
   !
END SUBROUTINE
!
SUBROUTINE test_cusolver_cmplx(ndim,mat)
   !
   USE cudafor
   USE cusolverdn
   !
   IMPLICIT NONE
   !
   INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14,200)
   !
   INTEGER, INTENT(IN) :: ndim
   COMPLEX(DP), INTENT(INOUT) :: mat(ndim,ndim)
   !
   INTEGER :: lwork
   INTEGER :: ierr
   INTEGER :: cr
   INTEGER :: t1
   INTEGER :: t2
   INTEGER, ALLOCATABLE :: ipiv(:)
   INTEGER, DEVICE :: info_d
   INTEGER, DEVICE, ALLOCATABLE :: ipiv_d(:)
   !
   COMPLEX(DP), ALLOCATABLE :: work(:)
   COMPLEX(DP), DEVICE, ALLOCATABLE :: mat_d(:,:)
   COMPLEX(DP), DEVICE, ALLOCATABLE :: work_d(:)
   !
   TYPE(cusolverDnHandle) :: cusolver_h
   !
   ierr = cusolverDnCreate(cusolver_h)
   ierr = cusolverDnZgetrf_bufferSize(cusolver_h,ndim,ndim,mat_d,ndim,lwork)
   !
   ALLOCATE(ipiv(ndim))
   ALLOCATE(work(1))
   ALLOCATE(mat_d(ndim,ndim))
   ALLOCATE(ipiv_d(ndim))
   ALLOCATE(work_d(lwork))
   !
   mat_d(:,:) = mat
   !
   ierr = cudaDeviceSynchronize()
   !
   CALL SYSTEM_CLOCK(COUNT_RATE=cr)
   CALL SYSTEM_CLOCK(t1)
   !
   ierr = cusolverDnZgetrf(cusolver_h,ndim,ndim,mat_d,ndim,work_d,ipiv_d,info_d)
   !
   CALL ZGETRI(ndim,mat,ndim,ipiv,work,-1,ierr)
   !
   lwork = CEILING(REAL(work(1)))
   !
   DEALLOCATE(work)
   ALLOCATE(work(lwork))
   !
   ipiv(:) = ipiv_d
   mat(:,:) = mat_d
   !
   CALL ZGETRI(ndim,mat,ndim,ipiv,work,lwork,ierr)
   !
   CALL SYSTEM_CLOCK(t2)
   !
   PRINT *,"time2:",REAL(t2-t1)/REAL(cr)
   !
   DEALLOCATE(ipiv)
   DEALLOCATE(work)
   DEALLOCATE(mat_d)
   DEALLOCATE(ipiv_d)
   DEALLOCATE(work_d)
   !
   ierr = cusolverDnDestroy(cusolver_h)
   !
END SUBROUTINE
