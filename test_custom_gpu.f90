SUBROUTINE test_custom_gpu_real(ndim,mat)
   !
   USE cudafor
   USE cublas_v2
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
   INTEGER :: lwork2
   INTEGER :: ierr
   INTEGER :: ii
   INTEGER :: jj
   INTEGER :: jp
   INTEGER :: cr
   INTEGER :: t1
   INTEGER :: t2
   INTEGER, DEVICE :: info_d
   INTEGER, ALLOCATABLE :: ipiv(:)
   INTEGER, DEVICE, ALLOCATABLE :: ipiv_d(:)
   !
   REAL(DP), DEVICE, ALLOCATABLE :: mat_d(:,:)
   REAL(DP), DEVICE, ALLOCATABLE :: work_d(:)
   !
   REAL(DP), PARAMETER :: zero = 0._DP
   REAL(DP), PARAMETER :: one = 1._DP
   !
   TYPE(cublasHandle) :: cublas_h
   TYPE(cusolverDnHandle) :: cusolver_h
   !
   ierr = cublasCreate(cublas_h)
   ierr = cusolverDnCreate(cusolver_h)
   ierr = cusolverDnDgetrf_bufferSize(cusolver_h,ndim,ndim,mat_d,ndim,lwork)
   ierr = cusolverDnDtrtri_buffersize(cusolver_h,CUBLAS_FILL_MODE_UPPER,&
   & CUBLAS_DIAG_NON_UNIT,ndim,mat_d,ndim,lwork2)
   !
   lwork = MAX(lwork,lwork2)
   lwork = MAX(lwork,ndim**2)
   !
   ALLOCATE(mat_d(ndim,ndim))
   ALLOCATE(ipiv_d(ndim))
   ALLOCATE(work_d(lwork))
   ALLOCATE(ipiv(ndim))
   !
   mat_d = mat
   !
   ierr = cudaDeviceSynchronize()
   !
   CALL SYSTEM_CLOCK(COUNT_RATE=cr)
   CALL SYSTEM_CLOCK(t1)
   !
   ierr = cusolverDnDgetrf(cusolver_h,ndim,ndim,mat_d,ndim,work_d,ipiv_d,info_d)
   !
   ipiv = ipiv_d
   !
   ierr = cusolverDnDtrtri(cusolver_h,CUBLAS_FILL_MODE_UPPER,&
   & CUBLAS_DIAG_NON_UNIT,ndim,mat_d,ndim,work_d,lwork,info_d)
   !
   DO jj = 1,ndim
      !$cuf kernel do(1) <<<*,*>>>
      DO ii = jj+1,ndim
         work_d(ii+(jj-1)*ndim) = mat_d(ii,jj)
         mat_d(ii,jj) = zero
      ENDDO
   ENDDO
   !
   ierr = cublasDtrsm_v2(cublas_h,CUBLAS_SIDE_RIGHT,CUBLAS_FILL_MODE_LOWER,&
   & CUBLAS_OP_N,CUBLAS_DIAG_UNIT,ndim,ndim,one,work_d,ndim,mat_d,ndim)
   !
   DO jj = ndim-1,1,-1
      jp = ipiv(jj)
      !
      IF(jp /= jj) THEN
         ierr = cublasDswap_v2(cublas_h,ndim,mat_d(1,jj),1,mat_d(1,jp),1)
      ENDIF
   ENDDO
   !
   ierr = cudaDeviceSynchronize()
   !
   CALL SYSTEM_CLOCK(t2)
   !
   PRINT *,"time5:",REAL(t2-t1)/REAL(cr)
   !
   mat = mat_d
   !
   DEALLOCATE(mat_d)
   DEALLOCATE(ipiv_d)
   DEALLOCATE(work_d)
   DEALLOCATE(ipiv)
   !
   ierr = cublasDestroy(cublas_h)
   ierr = cusolverDnDestroy(cusolver_h)
   !
END SUBROUTINE
!
SUBROUTINE test_custom_gpu_cmplx(ndim,mat)
   !
   USE cudafor
   USE cublas_v2
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
   INTEGER :: lwork2
   INTEGER :: ierr
   INTEGER :: ii
   INTEGER :: jj
   INTEGER :: jp
   INTEGER :: cr
   INTEGER :: t1
   INTEGER :: t2
   INTEGER, DEVICE :: info_d
   INTEGER, ALLOCATABLE :: ipiv(:)
   INTEGER, DEVICE, ALLOCATABLE :: ipiv_d(:)
   !
   COMPLEX(DP), DEVICE, ALLOCATABLE :: mat_d(:,:)
   COMPLEX(DP), DEVICE, ALLOCATABLE :: work_d(:)
   !
   COMPLEX(DP), PARAMETER :: zero = (0._DP,0._DP)
   COMPLEX(DP), PARAMETER :: one = (1._DP,0._DP)
   !
   TYPE(cublasHandle) :: cublas_h
   TYPE(cusolverDnHandle) :: cusolver_h
   !
   ierr = cublasCreate(cublas_h)
   ierr = cusolverDnCreate(cusolver_h)
   ierr = cusolverDnZgetrf_bufferSize(cusolver_h,ndim,ndim,mat_d,ndim,lwork)
   ierr = cusolverDnZtrtri_buffersize(cusolver_h,CUBLAS_FILL_MODE_UPPER,&
   & CUBLAS_DIAG_NON_UNIT,ndim,mat_d,ndim,lwork2)
   !
   lwork = MAX(lwork,lwork2)
   lwork = MAX(lwork,ndim**2)
   !
   ALLOCATE(mat_d(ndim,ndim))
   ALLOCATE(ipiv_d(ndim))
   ALLOCATE(work_d(lwork))
   ALLOCATE(ipiv(ndim))
   !
   mat_d = mat
   !
   ierr = cudaDeviceSynchronize()
   !
   CALL SYSTEM_CLOCK(COUNT_RATE=cr)
   CALL SYSTEM_CLOCK(t1)
   !
   ierr = cusolverDnZgetrf(cusolver_h,ndim,ndim,mat_d,ndim,work_d,ipiv_d,info_d)
   !
   ipiv = ipiv_d
   !
   ierr = cusolverDnZtrtri(cusolver_h,CUBLAS_FILL_MODE_UPPER,&
   & CUBLAS_DIAG_NON_UNIT,ndim,mat_d,ndim,work_d,lwork,info_d)
   !
   DO jj = 1,ndim
      !$cuf kernel do(1) <<<*,*>>>
      DO ii = jj+1,ndim
         work_d(ii+(jj-1)*ndim) = mat_d(ii,jj)
         mat_d(ii,jj) = zero
      ENDDO
   ENDDO
   !
   ierr = cublasZtrsm_v2(cublas_h,CUBLAS_SIDE_RIGHT,CUBLAS_FILL_MODE_LOWER,&
   & CUBLAS_OP_N,CUBLAS_DIAG_UNIT,ndim,ndim,one,work_d,ndim,mat_d,ndim)
   !
   DO jj = ndim-1,1,-1
      jp = ipiv(jj)
      !
      IF(jp /= jj) THEN
         ierr = cublasZswap_v2(cublas_h,ndim,mat_d(1,jj),1,mat_d(1,jp),1)
      ENDIF
   ENDDO
   !
   ierr = cudaDeviceSynchronize()
   !
   CALL SYSTEM_CLOCK(t2)
   !
   PRINT *,"time5:",REAL(t2-t1)/REAL(cr)
   !
   mat = mat_d
   !
   DEALLOCATE(mat_d)
   DEALLOCATE(ipiv_d)
   DEALLOCATE(work_d)
   DEALLOCATE(ipiv)
   !
   ierr = cublasDestroy(cublas_h)
   ierr = cusolverDnDestroy(cusolver_h)
   !
END SUBROUTINE
