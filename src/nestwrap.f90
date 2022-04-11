MODULE nestwrapper

!  Nested sampling includes
   USE params
   USE like
   IMPLICIT NONE

CONTAINS

! Wrapper around Likelihood Function
! Cube(1:n_dim) has nonphysical parameters
! scale Cube(1:n_dim) & return the scaled parameters in Cube(1:n_dim) &
! additional parameters in Cube(n_dim+1:nPar)
! return the log-likelihood in lnew

   SUBROUTINE getLogLike(Cube, n_dim, nPar, lnew, context)

      INTEGER                        :: n_dim, nPar, context
      REAL*8                         ::  lnew, Cube(nPar)
      INTEGER                        ::  i

      CALL FUserbuild(1.d0, i, lnew, i, i, i, Cube, i)
      IF (lnew <= -1.d10) THEN
         lnew = -1d90
      END IF

   END SUBROUTINE getLogLike

!-----*-----------------------------------------------------------------

! dumper, called after every updInt*10 iterations

   SUBROUTINE dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr, maxLogLike, logZ, logZerr, context)

      IMPLICIT NONE

      INTEGER                                 :: nSamples                ! number of samples in posterior array
      INTEGER                                 :: nlive                ! number of live points
      INTEGER                                 :: nPar                        ! number of parameters saved (physical plus derived)
      DOUBLE PRECISION, POINTER               :: physLive(:, :)        ! array containing the last set of live points
      DOUBLE PRECISION, POINTER               :: posterior(:, :)        ! array with the posterior distribution
      DOUBLE PRECISION, POINTER               :: paramConstr(:)        ! array with mean, sigmas, maxlike & MAP parameters
      DOUBLE PRECISION                        :: maxLogLike        ! max loglikelihood value
      DOUBLE PRECISION                        :: logZ                ! log evidence
      DOUBLE PRECISION                        :: logZerr                ! error on log evidence
      INTEGER                                 :: context         ! not required by MultiNest, any additional information user wants to pass

   END SUBROUTINE dumper

!-----*-----------------------------------------------------------------

END MODULE nestwrapper
