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
   INTEGER, DEVICE, ALLOCATABLE :: piv_d(:)
   !
   REAL(DP) :: tmp
   REAL(DP), DEVICE, ALLOCATABLE :: mat_d(:,:)
   REAL(DP), DEVICE, ALLOCATABLE :: work_d(:)
   !
   REAL(DP), PARAMETER :: zero = 0._DP
   REAL(DP), PARAMETER :: one = 1._DP
   !
   TYPE(cusolverDnHandle) :: cusolver_h
   !
   ierr = cusolverDnCreate(cusolver_h)
   ierr = cusolverDnDgetrf_bufferSize(cusolver_h,ndim,ndim,mat_d,ndim,lwork)
   ierr = cusolverDnDtrtri_buffersize(cusolver_h,CUBLAS_FILL_MODE_UPPER,&
   & CUBLAS_DIAG_NON_UNIT,ndim,mat_d,ndim,lwork2)
   !
   lwork = MAX(lwork,lwork2)
   lwork = MAX(lwork,ndim**2)
   !
   ALLOCATE(mat_d(ndim,ndim))
   ALLOCATE(piv_d(ndim))
   ALLOCATE(work_d(lwork))
   !
   mat_d(:,:) = mat
   !
   ierr = cudaDeviceSynchronize()
   !
   CALL SYSTEM_CLOCK(COUNT_RATE=cr)
   CALL SYSTEM_CLOCK(t1)
   !
   ierr = cusolverDnDgetrf(cusolver_h,ndim,ndim,mat_d,ndim,work_d,piv_d,info_d)
   ierr = cusolverDnDtrtri(cusolver_h,CUBLAS_FILL_MODE_UPPER,&
   & CUBLAS_DIAG_NON_UNIT,ndim,mat_d,ndim,work_d,lwork,info_d)
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
   DEALLOCATE(mat_d)
   DEALLOCATE(piv_d)
   DEALLOCATE(work_d)
   !
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
   INTEGER, DEVICE, ALLOCATABLE :: piv_d(:)
   !
   COMPLEX(DP) :: tmp
   COMPLEX(DP), DEVICE, ALLOCATABLE :: mat_d(:,:)
   COMPLEX(DP), DEVICE, ALLOCATABLE :: work_d(:)
   !
   COMPLEX(DP), PARAMETER :: zero = (0._DP,0._DP)
   COMPLEX(DP), PARAMETER :: one = (1._DP,0._DP)
   !
   TYPE(cusolverDnHandle) :: cusolver_h
   !
   ierr = cusolverDnCreate(cusolver_h)
   ierr = cusolverDnZgetrf_bufferSize(cusolver_h,ndim,ndim,mat_d,ndim,lwork)
   ierr = cusolverDnZtrtri_buffersize(cusolver_h,CUBLAS_FILL_MODE_UPPER,&
   & CUBLAS_DIAG_NON_UNIT,ndim,mat_d,ndim,lwork2)
   !
   lwork = MAX(lwork,lwork2)
   lwork = MAX(lwork,ndim**2)
   !
   ALLOCATE(mat_d(ndim,ndim))
   ALLOCATE(piv_d(ndim))
   ALLOCATE(work_d(lwork))
   !
   mat_d(:,:) = mat
   !
   ierr = cudaDeviceSynchronize()
   !
   CALL SYSTEM_CLOCK(COUNT_RATE=cr)
   CALL SYSTEM_CLOCK(t1)
   !
   ierr = cusolverDnZgetrf(cusolver_h,ndim,ndim,mat_d,ndim,work_d,piv_d,info_d)
   ierr = cusolverDnZtrtri(cusolver_h,CUBLAS_FILL_MODE_UPPER,&
   & CUBLAS_DIAG_NON_UNIT,ndim,mat_d,ndim,work_d,lwork,info_d)
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
   DEALLOCATE(mat_d)
   DEALLOCATE(piv_d)
   DEALLOCATE(work_d)
   !
   ierr = cusolverDnDestroy(cusolver_h)
   !
END SUBROUTINE
