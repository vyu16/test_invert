PROGRAM test_cusolver
   !
   USE cudafor
   USE cusolverdn
   !
   IMPLICIT NONE
   !
   INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14,200)
   !
   CHARACTER(100) :: arg1
   !
   INTEGER :: ndim
   INTEGER :: nrand
   INTEGER :: lwork
   INTEGER :: ierr
   INTEGER :: cr
   INTEGER :: t1
   INTEGER :: t2
   INTEGER :: ii
   INTEGER :: jj
   INTEGER, DEVICE :: info_d
   !
   INTEGER, ALLOCATABLE :: seed(:)
   INTEGER, ALLOCATABLE :: ipiv(:)
   INTEGER, DEVICE, ALLOCATABLE :: ipiv_d(:)
   !
   REAL(DP), ALLOCATABLE :: mat(:,:)
   REAL(DP), ALLOCATABLE :: work(:)
   REAL(DP), DEVICE, ALLOCATABLE :: mat_d(:,:)
   REAL(DP), DEVICE, ALLOCATABLE :: work_d(:)
   !
   TYPE(cusolverDnHandle) :: cusolver_h
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
   ierr = cusolverDnCreate(cusolver_h)
   ierr = cusolverDnDgetrf_bufferSize(cusolver_h,ndim,ndim,mat_d,ndim,lwork)
   !
   CALL RANDOM_SEED(SIZE=nrand)
   !
   ALLOCATE(seed(nrand))
   ALLOCATE(mat(ndim,ndim))
   ALLOCATE(ipiv(ndim))
   ALLOCATE(work(1))
   ALLOCATE(mat_d(ndim,ndim))
   ALLOCATE(ipiv_d(ndim))
   ALLOCATE(work_d(lwork))
   !
   seed = 123
   !
   CALL RANDOM_SEED(PUT=seed)
   CALL RANDOM_NUMBER(mat)
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
   CALL DGETRI(ndim,mat,ndim,ipiv,work,-1,ierr)
   !
   lwork = CEILING(work(1))
   !
   DEALLOCATE(seed)
   DEALLOCATE(work)
   ALLOCATE(work(lwork))
   !
   ipiv = ipiv_d
   mat = mat_d
   !
   CALL DGETRI(ndim,mat,ndim,ipiv,work,lwork,ierr)
   !
   CALL SYSTEM_CLOCK(t2)
   !
   DO jj = 1,3
      DO ii = 1,3
         WRITE(*,'(E18.8)') mat(ii,jj)
      ENDDO
   ENDDO
   !
   PRINT *,"Time:",REAL(t2-t1)/REAL(cr)
   !
   DEALLOCATE(mat)
   DEALLOCATE(ipiv)
   DEALLOCATE(work)
   DEALLOCATE(mat_d)
   DEALLOCATE(ipiv_d)
   DEALLOCATE(work_d)
   !
   ierr = cusolverDnDestroy(cusolver_h)
   !
END PROGRAM
