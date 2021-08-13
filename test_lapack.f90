PROGRAM test_cpu
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
   !
   INTEGER, ALLOCATABLE :: seed(:)
   INTEGER, ALLOCATABLE :: ipiv(:)
   !
   REAL(DP), ALLOCATABLE :: mat(:,:)
   REAL(DP), ALLOCATABLE :: work(:)
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
   ALLOCATE(ipiv(ndim))
   ALLOCATE(work(1))
   !
   seed = 123
   !
   CALL RANDOM_SEED(PUT=seed)
   CALL RANDOM_NUMBER(mat)
   !
   CALL SYSTEM_CLOCK(COUNT_RATE=cr)
   CALL SYSTEM_CLOCK(t1)
   !
   CALL DGETRF(ndim,ndim,mat,ndim,ipiv,ierr)
   !
   CALL DGETRI(ndim,mat,ndim,ipiv,work,-1,ierr)
   !
   lwork = CEILING(work(1))
   !
   DEALLOCATE(seed)
   DEALLOCATE(work)
   ALLOCATE(work(lwork))
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
   !
END PROGRAM
