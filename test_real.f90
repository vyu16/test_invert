PROGRAM test_invert_real
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
   REAL(DP), ALLOCATABLE :: mat1(:,:)
   REAL(DP), ALLOCATABLE :: mat2(:,:)
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
   ALLOCATE(mat1(ndim,ndim))
   ALLOCATE(mat2(ndim,ndim))
   !
   seed(:) = 123
   !
   CALL RANDOM_SEED(PUT=seed)
   CALL RANDOM_NUMBER(mat1)
   !
   mat2(:,:) = mat1
   !
   CALL test_cublas_real(ndim,mat1)
   CALL test_custom_gpu_real(ndim,mat2)
   !
   PRINT *,"error:",MAXVAL(ABS(mat2-mat1))
   !
   DEALLOCATE(seed)
   DEALLOCATE(mat1)
   DEALLOCATE(mat2)
   !
END PROGRAM
