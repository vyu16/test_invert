SUBROUTINE test_custom_real(ndim,mat)
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
   lwork = ndim*ndim
   !
   ALLOCATE(ipiv(ndim))
   ALLOCATE(work(lwork))
   !
   CALL SYSTEM_CLOCK(COUNT_RATE=cr)
   CALL SYSTEM_CLOCK(t1)
   !
   CALL DGETRF(ndim,ndim,mat,ndim,ipiv,ierr)
   CALL DGETRI2(ndim,mat,ndim,ipiv,work,lwork,ierr)
   !
   CALL SYSTEM_CLOCK(t2)
   !
   PRINT *,"time4:",REAL(t2-t1)/REAL(cr)
   !
   DEALLOCATE(ipiv)
   DEALLOCATE(work)
   !
END SUBROUTINE
!
SUBROUTINE test_custom_cmplx(ndim,mat)
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
   lwork = ndim*ndim
   !
   ALLOCATE(ipiv(ndim))
   ALLOCATE(work(lwork))
   !
   CALL SYSTEM_CLOCK(COUNT_RATE=cr)
   CALL SYSTEM_CLOCK(t1)
   !
   CALL ZGETRF(ndim,ndim,mat,ndim,ipiv,ierr)
   CALL ZGETRI2(ndim,mat,ndim,ipiv,work,lwork,ierr)
   !
   CALL SYSTEM_CLOCK(t2)
   !
   PRINT *,"time4:",REAL(t2-t1)/REAL(cr)
   !
   DEALLOCATE(ipiv)
   DEALLOCATE(work)
   !
END SUBROUTINE
!
SUBROUTINE DGETRI2(ndim,mat,lda,ipiv,work,lwork,info)
   !
   IMPLICIT NONE
   !
   INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14,200)
   !
   INTEGER, INTENT(IN) :: ndim
   INTEGER, INTENT(IN) :: lda
   INTEGER, INTENT(IN) :: ipiv(ndim)
   INTEGER, INTENT(IN) :: lwork
   INTEGER, INTENT(OUT) :: info
   REAL(DP), INTENT(INOUT) :: mat(ndim,ndim)
   REAL(DP), INTENT(OUT) :: work(ndim*ndim)
   !
   INTEGER :: ii
   INTEGER :: jj
   INTEGER :: jp
   !
   REAL(DP), PARAMETER :: zero = 0._DP
   REAL(DP), PARAMETER :: one = 1._DP
   !
   CALL DTRTRI('U','N',ndim,mat,lda,info)
   !
   DO jj = 1,ndim
      DO ii = jj+1,ndim
         work(ii+(jj-1)*ndim) = mat(ii,jj)
         mat(ii,jj) = zero
      ENDDO
   ENDDO
   !
   CALL DTRSM('R','L','N','U',ndim,ndim,one,work,ndim,mat,lda)
   !
   DO jj = ndim-1,1,-1
      jp = ipiv(jj)
      !
      IF(jp /= jj) THEN
         CALL DSWAP(ndim,mat(1,jj),1,mat(1,jp),1)
      ENDIF
   ENDDO
   !
END SUBROUTINE
!
SUBROUTINE ZGETRI2(ndim,mat,lda,ipiv,work,lwork,info)
   !
   IMPLICIT NONE
   !
   INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14,200)
   !
   INTEGER, INTENT(IN) :: ndim
   INTEGER, INTENT(IN) :: lda
   INTEGER, INTENT(IN) :: ipiv(ndim)
   INTEGER, INTENT(IN) :: lwork
   INTEGER, INTENT(OUT) :: info
   COMPLEX(DP), INTENT(INOUT) :: mat(ndim,ndim)
   COMPLEX(DP), INTENT(OUT) :: work(ndim*ndim)
   !
   INTEGER :: ii
   INTEGER :: jj
   INTEGER :: jp
   !
   COMPLEX(DP), PARAMETER :: zero = (0._DP,0._DP)
   COMPLEX(DP), PARAMETER :: one = (1._DP,0._DP)
   !
   CALL ZTRTRI('U','N',ndim,mat,lda,info)
   !
   DO jj = ndim,1,-1
      DO ii = jj+1,ndim
         work(ii+(jj-1)*ndim) = mat(ii,jj)
         mat(ii,jj) = zero
      ENDDO
   ENDDO
   !
   CALL ZTRSM('R','L','N','U',ndim,ndim,one,work,ndim,mat,lda)
   !
   DO jj = ndim-1,1,-1
      jp = ipiv(jj)
      !
      IF(jp /= jj) THEN
         CALL ZSWAP(ndim,mat(1,jj),1,mat(1,jp),1)
      ENDIF
   ENDDO
   !
END SUBROUTINE
