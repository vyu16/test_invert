SUBROUTINE test_custom_gpu_real(ndim,mat)
   !
   USE cudafor
   USE cublas
   USE cusolverdn
   !
   IMPLICIT NONE
   !
   INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14,200)
   !
   INTEGER, INTENT(IN) :: ndim
   REAL(DP), INTENT(INOUT) :: mat(ndim,ndim)
   !
   INTEGER(8) :: ndim_i8
   INTEGER(8) :: lwork
   INTEGER(8) :: lwork_d
   INTEGER(8) :: lwork2
   INTEGER(8) :: lwork2_d
   INTEGER :: ierr
   INTEGER :: ii
   INTEGER :: jj
   INTEGER :: jp
   INTEGER :: cr
   INTEGER :: t1
   INTEGER :: t2
   INTEGER, DEVICE :: info_d
   INTEGER(8), DEVICE, ALLOCATABLE :: piv_d(:)
   REAL(DP) :: tmp
   REAL(DP), ALLOCATABLE :: work(:)
   REAL(DP), DEVICE, ALLOCATABLE :: work_d(:)
   REAL(DP), DEVICE, ALLOCATABLE :: mat_d(:,:)
   !
   REAL(DP), PARAMETER :: zero = 0._DP
   REAL(DP), PARAMETER :: one = 1._DP
   !
   TYPE(cusolverDnHandle) :: cusolver_h
   TYPE(cusolverDnParams) :: cusolver_p
   !
   ndim_i8 = ndim
   !
   ierr = cusolverDnCreate(cusolver_h)
   ierr = cusolverDnCreateParams(cusolver_p)
   !
   ierr = cusolverDnXgetrf_bufferSize(cusolver_h,cusolver_p,ndim_i8,ndim_i8,&
        & cudaDataType(CUDA_R_64F),mat_d,ndim_i8,cudaDataType(CUDA_R_64F),&
        & lwork_d,lwork)
   ierr = cusolverDnXtrtri_buffersize(cusolver_h,CUBLAS_FILL_MODE_UPPER,&
        & CUBLAS_DIAG_NON_UNIT,ndim_i8,cudaDataType(CUDA_R_64F),mat_d,ndim_i8,&
        & lwork2_d,lwork2)
   !
   lwork = MAX(lwork,lwork2,8*ndim_i8**2)
   lwork_d = MAX(lwork_d,lwork2_d,8*ndim_i8**2)
   !
   ALLOCATE(work(lwork/8))
   ALLOCATE(work_d(lwork_d/8))
   ALLOCATE(mat_d(ndim,ndim))
   ALLOCATE(piv_d(ndim))
   !
   mat_d(:,:) = mat
   !
   ierr = cudaDeviceSynchronize()
   !
   CALL SYSTEM_CLOCK(COUNT_RATE=cr)
   CALL SYSTEM_CLOCK(t1)
   !
   ierr = cusolverDnXgetrf(cusolver_h,cusolver_p,ndim_i8,ndim_i8,&
        & cudaDataType(CUDA_R_64F),mat_d,ndim_i8,piv_d,&
        & cudaDataType(CUDA_R_64F),work_d,lwork_d,work,lwork,info_d)
   ierr = cusolverDnXtrtri(cusolver_h,CUBLAS_FILL_MODE_UPPER,&
        & CUBLAS_DIAG_NON_UNIT,ndim_i8,cudaDataType(CUDA_R_64F),mat_d,ndim_i8,&
        & work_d,lwork_d,work,lwork,info_d)
   !
   !$acc parallel loop collapse(2)
   DO jj = 1,ndim-1
      DO ii = 2,ndim
         IF(ii > jj) THEN
            work_d(ii+(jj-1)*ndim) = mat_d(ii,jj)
            mat_d(ii,jj) = zero
         ENDIF
      ENDDO
   ENDDO
   !$acc end parallel
   !
   CALL DTRSM('R','L','N','U',ndim,ndim,one,work_d,ndim,mat_d,ndim)
   !
   !$acc parallel
   !$acc loop seq
   DO jj = ndim-1,1,-1
      jp = piv_d(jj)
      !
      IF(jp /= jj) THEN
         !$acc loop
         DO ii = 1,ndim
            tmp = mat_d(ii,jj)
            mat_d(ii,jj) = mat_d(ii,jp)
            mat_d(ii,jp) = tmp
         ENDDO
      ENDIF
   ENDDO
   !$acc end parallel
   !
   ierr = cudaDeviceSynchronize()
   !
   CALL SYSTEM_CLOCK(t2)
   !
   PRINT *,"time5:",REAL(t2-t1)/REAL(cr)
   !
   mat(:,:) = mat_d
   !
   DEALLOCATE(work)
   DEALLOCATE(work_d)
   DEALLOCATE(mat_d)
   DEALLOCATE(piv_d)
   !
   ierr = cusolverDnDestroyParams(cusolver_p)
   ierr = cusolverDnDestroy(cusolver_h)
   !
END SUBROUTINE
!
SUBROUTINE test_custom_gpu_cmplx(ndim,mat)
   !
   USE cudafor
   USE cublas
   USE cusolverdn
   !
   IMPLICIT NONE
   !
   INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14,200)
   !
   INTEGER, INTENT(IN) :: ndim
   COMPLEX(DP), INTENT(INOUT) :: mat(ndim,ndim)
   !
   INTEGER(8) :: ndim_i8
   INTEGER(8) :: lwork
   INTEGER(8) :: lwork_d
   INTEGER(8) :: lwork2
   INTEGER(8) :: lwork2_d
   INTEGER :: ierr
   INTEGER :: ii
   INTEGER :: jj
   INTEGER :: jp
   INTEGER :: cr
   INTEGER :: t1
   INTEGER :: t2
   INTEGER, DEVICE :: info_d
   INTEGER(8), DEVICE, ALLOCATABLE :: piv_d(:)
   COMPLEX(DP) :: tmp
   COMPLEX(DP), ALLOCATABLE :: work(:)
   COMPLEX(DP), DEVICE, ALLOCATABLE :: work_d(:)
   COMPLEX(DP), DEVICE, ALLOCATABLE :: mat_d(:,:)
   !
   COMPLEX(DP), PARAMETER :: zero = (0._DP,0._DP)
   COMPLEX(DP), PARAMETER :: one = (1._DP,0._DP)
   !
   TYPE(cusolverDnHandle) :: cusolver_h
   TYPE(cusolverDnParams) :: cusolver_p
   !
   ndim_i8 = ndim
   !
   ierr = cusolverDnCreate(cusolver_h)
   ierr = cusolverDnCreateParams(cusolver_p)
   !
   ierr = cusolverDnXgetrf_bufferSize(cusolver_h,cusolver_p,ndim_i8,ndim_i8,&
        & cudaDataType(CUDA_C_64F),mat_d,ndim_i8,cudaDataType(CUDA_C_64F),&
        & lwork_d,lwork)
   ierr = cusolverDnXtrtri_buffersize(cusolver_h,CUBLAS_FILL_MODE_UPPER,&
        & CUBLAS_DIAG_NON_UNIT,ndim_i8,cudaDataType(CUDA_C_64F),mat_d,ndim_i8,&
        & lwork2_d,lwork2)
   !
   lwork = MAX(lwork,lwork2,16*ndim_i8**2)
   lwork_d = MAX(lwork_d,lwork2_d,16*ndim_i8**2)
   !
   ALLOCATE(work(lwork/16))
   ALLOCATE(work_d(lwork_d/16))
   ALLOCATE(mat_d(ndim,ndim))
   ALLOCATE(piv_d(ndim))
   !
   mat_d(:,:) = mat
   !
   ierr = cudaDeviceSynchronize()
   !
   CALL SYSTEM_CLOCK(COUNT_RATE=cr)
   CALL SYSTEM_CLOCK(t1)
   !
   ierr = cusolverDnXgetrf(cusolver_h,cusolver_p,ndim_i8,ndim_i8,&
        & cudaDataType(CUDA_C_64F),mat_d,ndim_i8,piv_d,&
        & cudaDataType(CUDA_C_64F),work_d,lwork_d,work,lwork,info_d)
   ierr = cusolverDnXtrtri(cusolver_h,CUBLAS_FILL_MODE_UPPER,&
        & CUBLAS_DIAG_NON_UNIT,ndim_i8,cudaDataType(CUDA_C_64F),mat_d,ndim_i8,&
        & work_d,lwork_d,work,lwork,info_d)
   !
   !$acc parallel loop collapse(2)
   DO jj = 1,ndim-1
      DO ii = 2,ndim
         IF(ii > jj) THEN
            work_d(ii+(jj-1)*ndim) = mat_d(ii,jj)
            mat_d(ii,jj) = zero
         ENDIF
      ENDDO
   ENDDO
   !$acc end parallel
   !
   CALL ZTRSM('R','L','N','U',ndim,ndim,one,work_d,ndim,mat_d,ndim)
   !
   !$acc parallel
   !$acc loop seq
   DO jj = ndim-1,1,-1
      jp = piv_d(jj)
      !
      IF(jp /= jj) THEN
         !$acc loop
         DO ii = 1,ndim
            tmp = mat_d(ii,jj)
            mat_d(ii,jj) = mat_d(ii,jp)
            mat_d(ii,jp) = tmp
         ENDDO
      ENDIF
   ENDDO
   !$acc end parallel
   !
   ierr = cudaDeviceSynchronize()
   !
   CALL SYSTEM_CLOCK(t2)
   !
   PRINT *,"time5:",REAL(t2-t1)/REAL(cr)
   !
   mat(:,:) = mat_d
   !
   DEALLOCATE(work)
   DEALLOCATE(work_d)
   DEALLOCATE(mat_d)
   DEALLOCATE(piv_d)
   !
   ierr = cusolverDnDestroyParams(cusolver_p)
   ierr = cusolverDnDestroy(cusolver_h)
   !
END SUBROUTINE
