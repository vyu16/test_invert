SUBROUTINE test_cublas_real(ndim,mat)
   !
   USE cudafor
   USE cublas_v2
   !
   IMPLICIT NONE
   !
   INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14,200)
   !
   INTEGER, INTENT(IN) :: ndim
   REAL(DP), INTENT(INOUT) :: mat(ndim,ndim)
   !
   INTEGER :: ierr
   INTEGER :: cr
   INTEGER :: t1
   INTEGER :: t2
   INTEGER, DEVICE :: info_d(1)
   INTEGER, DEVICE, ALLOCATABLE :: ipiv_d(:)
   !
   REAL(DP), DEVICE, ALLOCATABLE :: mat_d(:,:)
   REAL(DP), DEVICE, ALLOCATABLE :: work_d(:,:)
   !
   TYPE(C_DEVPTR), DEVICE :: mat_d_ptr(1)
   TYPE(C_DEVPTR), DEVICE :: work_d_ptr(1)
   !
   TYPE(cublasHandle) :: cublas_h
   !
   ALLOCATE(mat_d(ndim,ndim))
   ALLOCATE(ipiv_d(ndim))
   ALLOCATE(work_d(ndim,ndim))
   !
   mat_d = mat
   !
   ierr = cublasCreate(cublas_h)
   ierr = cudaDeviceSynchronize()
   !
   CALL SYSTEM_CLOCK(COUNT_RATE=cr)
   CALL SYSTEM_CLOCK(t1)
   !
   mat_d_ptr(1) = C_DEVLOC(mat_d(1,1))
   work_d_ptr(1) = C_DEVLOC(work_d(1,1))
   !
   ierr = cublasDgetrfBatched(cublas_h,ndim,mat_d_ptr,ndim,ipiv_d,info_d,1)
   ierr = cublasDgetriBatched(cublas_h,ndim,mat_d_ptr,ndim,ipiv_d,work_d_ptr,ndim,info_d,1)
   ierr = cudaDeviceSynchronize()
   !
   CALL SYSTEM_CLOCK(t2)
   !
   mat = work_d
   !
   PRINT *,"time3:",REAL(t2-t1)/REAL(cr)
   !
   DEALLOCATE(mat_d)
   DEALLOCATE(ipiv_d)
   DEALLOCATE(work_d)
   !
   ierr = cublasDestroy(cublas_h)
   !
END SUBROUTINE
!
SUBROUTINE test_cublas_cmplx(ndim,mat)
   !
   USE cudafor
   USE cublas_v2
   !
   IMPLICIT NONE
   !
   INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14,200)
   !
   INTEGER, INTENT(IN) :: ndim
   COMPLEX(DP), INTENT(INOUT) :: mat(ndim,ndim)
   !
   INTEGER :: ierr
   INTEGER :: cr
   INTEGER :: t1
   INTEGER :: t2
   INTEGER, DEVICE :: info_d(1)
   INTEGER, DEVICE, ALLOCATABLE :: ipiv_d(:)
   !
   COMPLEX(DP), DEVICE, ALLOCATABLE :: mat_d(:,:)
   COMPLEX(DP), DEVICE, ALLOCATABLE :: work_d(:,:)
   !
   TYPE(C_DEVPTR), DEVICE :: mat_d_ptr(1)
   TYPE(C_DEVPTR), DEVICE :: work_d_ptr(1)
   !
   TYPE(cublasHandle) :: cublas_h
   !
   ALLOCATE(mat_d(ndim,ndim))
   ALLOCATE(ipiv_d(ndim))
   ALLOCATE(work_d(ndim,ndim))
   !
   mat_d = mat
   !
   ierr = cublasCreate(cublas_h)
   ierr = cudaDeviceSynchronize()
   !
   CALL SYSTEM_CLOCK(COUNT_RATE=cr)
   CALL SYSTEM_CLOCK(t1)
   !
   mat_d_ptr(1) = C_DEVLOC(mat_d(1,1))
   work_d_ptr(1) = C_DEVLOC(work_d(1,1))
   !
   ierr = cublasZgetrfBatched(cublas_h,ndim,mat_d_ptr,ndim,ipiv_d,info_d,1)
   ierr = cublasZgetriBatched(cublas_h,ndim,mat_d_ptr,ndim,ipiv_d,work_d_ptr,ndim,info_d,1)
   ierr = cudaDeviceSynchronize()
   !
   CALL SYSTEM_CLOCK(t2)
   !
   mat = work_d
   !
   PRINT *,"time3:",REAL(t2-t1)/REAL(cr)
   !
   DEALLOCATE(mat_d)
   DEALLOCATE(ipiv_d)
   DEALLOCATE(work_d)
   !
   ierr = cublasDestroy(cublas_h)
   !
END SUBROUTINE
