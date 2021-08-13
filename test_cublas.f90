PROGRAM test_cublas
   !
   USE cudafor
   USE cublas_v2
   !
   IMPLICIT NONE
   !
   INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14,200)
   !
   CHARACTER(100) :: arg1
   !
   INTEGER :: ndim
   INTEGER :: nrand
   INTEGER :: ierr
   INTEGER :: cr
   INTEGER :: t1
   INTEGER :: t2
   INTEGER :: ii
   INTEGER :: jj
   INTEGER, DEVICE :: info_d(1)
   !
   INTEGER, ALLOCATABLE :: seed(:)
   INTEGER, DEVICE, ALLOCATABLE :: ipiv_d(:)
   !
   REAL(DP), ALLOCATABLE :: mat(:,:)
   REAL(DP), DEVICE, ALLOCATABLE :: mat_d(:,:)
   REAL(DP), DEVICE, ALLOCATABLE :: work_d(:,:)
   !
   TYPE(C_DEVPTR), DEVICE :: mat_d_ptr(1)
   TYPE(C_DEVPTR), DEVICE :: work_d_ptr(1)
   TYPE(cublasHandle) :: cublas_h
   !
   ndim = -1
   !
   IF(COMMAND_ARGUMENT_COUNT() == 1) THEN
      CALL GET_COMMAND_ARGUMENT(1,arg1)
      !
      READ(arg1,*) ndim
   ENDIF
   !
   ndim = MAX(ndim,10)
   !
   CALL RANDOM_SEED(SIZE=nrand)
   !
   ALLOCATE(seed(nrand))
   ALLOCATE(mat(ndim,ndim))
   ALLOCATE(mat_d(ndim,ndim))
   ALLOCATE(ipiv_d(ndim))
   ALLOCATE(work_d(ndim,ndim))
   !
   seed = 123
   !
   CALL RANDOM_SEED(PUT=seed)
   CALL RANDOM_NUMBER(mat)
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
   DO jj = 1,3
      DO ii = 1,3
         WRITE(*,'(E18.8)') mat(ii,jj)
      ENDDO
   ENDDO
   !
   PRINT *,"Time:",REAL(t2-t1)/REAL(cr)
   !
   DEALLOCATE(seed)
   DEALLOCATE(mat)
   DEALLOCATE(mat_d)
   DEALLOCATE(ipiv_d)
   DEALLOCATE(work_d)
   !
   ierr = cublasDestroy(cublas_h)
   !
END PROGRAM
