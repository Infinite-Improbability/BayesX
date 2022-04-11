MODULE matrix_utils

CONTAINS

   SUBROUTINE Transpose(nx, ny, M, T)

      INTEGER                        ::  nx, ny, i, j
      REAL*8                         ::  M(nx, ny), T(ny, nx)

! T_ij = M_ji

      DO j = 1, ny
         DO i = 1, nx
            T(j, i) = M(i, j)
         END DO
      END DO

      RETURN
   END SUBROUTINE Transpose
!=======================================================================

   SUBROUTINE MatrixProduct(nx1, ny1, M1, nx2, ny2, M2, M3)

      INTEGER                        ::  nx1, ny1, nx2, ny2, i, j, k
      REAL*8                         :: M1(nx1, ny1), M2(nx2, ny2), M3(nx1, ny2), sum

      IF (ny1 .ne. nx2) then
         if (myID == 0) write (*, *) 'ERROR in MatrixProduct: ny1 != nx2'
         stop
      end if

! M3_ij = M1_ik * M2_kj

      DO j = 1, ny2
         DO i = 1, nx1

            sum = 0.0
            DO k = 1, ny1
               IF (M1(i, k) .ne. 0.0 .and. M2(k, j) .ne. 0.0) then
                  sum = sum + M1(i, k)*M2(k, j)
               END IF
            END DO
            M3(i, j) = sum

         END DO
      END DO

      RETURN
   END SUBROUTINE MatrixProduct

!=======================================================================

END MODULE matrix_utils
