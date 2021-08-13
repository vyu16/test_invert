SUBROUTINE test_lapack_real(ndim,mat)
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
   !
   REAL(DP), ALLOCATABLE :: work(:)
   !
   ALLOCATE(ipiv(ndim))
   ALLOCATE(work(1))
   !
   CALL SYSTEM_CLOCK(COUNT_RATE=cr)
   CALL SYSTEM_CLOCK(t1)
   !
   CALL DGETRF(ndim,ndim,mat,ndim,ipiv,ierr)
   CALL DGETRI(ndim,mat,ndim,ipiv,work,-1,ierr)
   !
   lwork = CEILING(work(1))
   !
   DEALLOCATE(work)
   ALLOCATE(work(lwork))
   !
   CALL DGETRI(ndim,mat,ndim,ipiv,work,lwork,ierr)
   !
   CALL SYSTEM_CLOCK(t2)
   !
   PRINT *,"time1:",REAL(t2-t1)/REAL(cr)
   !
   DEALLOCATE(ipiv)
   DEALLOCATE(work)
   !
END SUBROUTINE
!
SUBROUTINE test_lapack_cmplx(ndim,mat)
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
   !
   COMPLEX(DP), ALLOCATABLE :: work(:)
   !
   ALLOCATE(ipiv(ndim))
   ALLOCATE(work(1))
   !
   CALL SYSTEM_CLOCK(COUNT_RATE=cr)
   CALL SYSTEM_CLOCK(t1)
   !
   CALL ZGETRF(ndim,ndim,mat,ndim,ipiv,ierr)
   CALL ZGETRI(ndim,mat,ndim,ipiv,work,-1,ierr)
   !
   lwork = CEILING(REAL(work(1)))
   !
   DEALLOCATE(work)
   ALLOCATE(work(lwork))
   !
   CALL ZGETRI(ndim,mat,ndim,ipiv,work,lwork,ierr)
   !
   CALL SYSTEM_CLOCK(t2)
   !
   PRINT *,"time1:",REAL(t2-t1)/REAL(cr)
   !
   DEALLOCATE(ipiv)
   DEALLOCATE(work)
   !
END SUBROUTINE
