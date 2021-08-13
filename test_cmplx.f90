PROGRAM test_invert_cmplx
   !
   IMPLICIT NONE
   !
   INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14,200)
   !
   CHARACTER(100) :: arg1
   !
   INTEGER :: ndim
   INTEGER :: nrand
   INTEGER, ALLOCATABLE :: seed(:)
   !
   REAL(DP), ALLOCATABLE :: temp(:,:)
   !
   COMPLEX(DP), ALLOCATABLE :: mat1(:,:)
   COMPLEX(DP), ALLOCATABLE :: mat2(:,:)
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
   ALLOCATE(temp(ndim,ndim))
   ALLOCATE(mat1(ndim,ndim))
   ALLOCATE(mat2(ndim,ndim))
   !
   seed(:) = 123
   !
   CALL RANDOM_SEED(PUT=seed)
   CALL RANDOM_NUMBER(temp)
   !
   mat1(:,:) = CMPLX(temp-1,temp+1,KIND=DP)
   mat2(:,:) = mat1
   !
   CALL test_lapack_cmplx(ndim,mat1)
!   CALL test_cusolver_cmplx(ndim,mat2)
   !
   PRINT *,"error:",MAXVAL(ABS(mat2-mat1))
   !
   DEALLOCATE(seed)
   DEALLOCATE(temp)
   DEALLOCATE(mat1)
   DEALLOCATE(mat2)
   !
END PROGRAM
