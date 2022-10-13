PROGRAM test_macro
   !
   IMPLICIT NONE
   !
#if defined(__PGI)
   PRINT *,"PGI"
#endif
   !
#if defined(__PGI) && (__PGIC__ > 22 || (__PGIC__ == 22 && __PGIC_MINOR__ > 7))
   PRINT *,"PGI >= 22.9"
#endif
   !
END PROGRAM
