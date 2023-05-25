MODULE like
   USE params
   USE constants
   USE CheckPars1
   USE ReadWrite
   USE MassModels
   USE GasModels
   USE priors1
   USE massfunction
   USE utilities
   USE RandomNS

   IMPLICIT NONE

CONTAINS

   SUBROUTINE FUserbuild(Cool, Nsystem, Lhood, Nstore, idummy, ndummy, Cube1, retcode)

      IMPLICIT NONE

#ifdef MPI
      INCLUDE 'mpif.h'
#endif

      INTEGER                         :: Nsystem, Nstore, idummy, ndummy, retcode, i, k, i1, j, flag, list(NAtoms), id, err
      REAL*8                          :: Cool, Lhood, Cube1(:), Cube(tot_dim), Cube2(tot_dim)
      REAL*8                          :: XRAYLhood, urv

!-----------------------------------------------------------------------
      flag = 0
      retcode = 1
      Cube(1:edim) = Cube1(1:edim)
      Cube(edim + 1:tot_dim) = 0d0
      XRAYLhood = 0.d0

!        add the parameters with delta priors

      DO i = 1, tot_dim
         CALL PClass(i, j)
         IF (j == 0) Cube(i + 1:tot_dim) = Cube(i:tot_dim - 1)
      END DO

!         Additional prior constraints: check parameters and return both an
!         additional contribution to the likelihood and a flag value (for
!         truncated distributions)

      CALL CheckPars(Cube, flag)
      IF (flag == 1) GOTO 999

#ifdef MPI
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, id, err)
#endif

!     Make predicted data:

      CALL PredictXrayCounts(Cube, flag)
      IF (flag == 1) GOTO 999

      ! Write predicted counts to file
      if (clusterDims == 0) THEN
         open (unit=105, form='formatted', file=trim(n_root)//'generated-data.txt', status='replace')
         do i = 1, size(xrayCpred)
            write (105, *) xrayCpred(i)
         end do
         close (105, status='KEEP')
      end if

!     Calculate Likelihood

      CALL XrayCountsLhood(XRAYLhood)
      Lhood = XRAYLhood

999   IF (flag == 0) THEN
         retcode = 1
      ELSE
         retcode = 0
         Lhood = B_A_D
         RETURN
      END IF

!        return the scaled (physical) parameters for nested sampling
      DO i = 1, tot_dim
         CALL Rescale_nest(Cube(i), i, i1)
      END DO

      list = 0
      list(NAtoms) = 1
      DO i = 2, NAtoms
         DO j = 1, NAtoms
            IF (list(j) == 0) CYCLE
            IF (Cube((i - 1)*NPars + 1) < Cube((list(j) - 1)*NPars + 1)) THEN
               k = j - 1
               exit
            END IF

            IF (j == NAtoms) THEN
               k = NAtoms
               exit
            END IF
         END DO

         list(1:k - 1) = list(2:k)
         list(k) = i
      END DO

      DO i = 1, NAtoms
         DO j = 1, NPars
            Cube2((i - 1)*NPars + j) = Cube((list(i) - 1)*NPars + j)
         END DO
      END DO
      Cube2(NAtoms*NPars + 1:tot_dim) = Cube(NAtoms*NPars + 1:tot_dim)
      Cube = Cube2

!        remove the parameters with delta priors

      k = 0
      DO i = 1, tot_dim
         CALL PClass(i, j)
         IF (j /= 0) THEN
            k = k + 1
            Cube1(k) = Cube(i)
         END IF
      END DO

      DO i = 1, NAtoms
         Cube1(k + 1:k + aux_dim) = aux(i, 1:aux_dim)
         k = k + aux_dim
      END DO

   END SUBROUTINE FUserbuild

!=======================================================================
   SUBROUTINE PredictXrayCounts(cube, flag)

      IMPLICIT NONE

      INTEGER                        ::  i, flag
      REAL*8                          ::  Cube(*)
!        initialize working arrays
      xrayCpred = 0.0d0
      xrayCmap = 0.d0
      DO i = 1, NAtoms
         CALL MakeXrayGasDistributions(i, flag)
         IF (flag == 1) THEN
            RETURN
         ELSEIF (flag == 2) THEN
            flag = 0
         END IF
      END DO

      RETURN

   END SUBROUTINE PredictXrayCounts

!=======================================================================
   SUBROUTINE XraycountsLhood(XRAYLhood)

      IMPLICIT NONE

      INTEGER                         ::  m
      REAL*8                           ::   XRAYLhood, sum

      sum = 0.0d0

      DO m = 1, LENx

         IF (xrayMask(m) .EQ. 0) THEN

            IF (xrayCpred(m) .EQ. 0d0) THEN
               sum = sum + (xrayCobs(m)*DLOG(xrayBG(m)) - xrayBG(m)) + (xrayBG_obs(m)*DLOG(bexpotime*xrayBG(m)/sexpotime) - (bexpotime*xrayBG(m)/sexpotime))
            ELSE
               sum = sum + (xrayCobs(m)*DLOG(xrayCpred(m)) - xrayCpred(m)) + (xrayBG_obs(m)*DLOG(bexpotime*xrayBG(m)/sexpotime) - (bexpotime*xrayBG(m)/sexpotime))
            END IF

         END IF

      END DO
      XRAYLhood = sum + XRAYLhood0
      ! WRITE(*,*)' XRAYLhood =', XRAYLhood

   END SUBROUTINE XraycountsLhood

!=======================================================================
   integer function Initialise(index)

      IMPLICIT NONE
      INTEGER                         ::  i, j, k, m, index, iend, idx(1), iostatus
      REAL*8                          :: Lhood, XRAYLhood, log_Cobs_factorial, log_BGobs_factorial, nullev, angfactor
      REAL*8                          :: Mmin, Mmax, d1, d2
      real*8                          :: M200_max
      CHARACTER(LEN=100)    :: string

      Initialise = 0
      tot_dim = NDim
      tot_atoms = NAtoms

!        reduce the dimensionality if there are parameters with delta priors
      edim = 0
      DO i = 1, tot_dim
         CALL PClass(i, j)
         IF (j > 0) edim = edim + 1
      END DO

!       set the total number of parameters

      !n_totPar = tot_dim + aux_dim*NAtoms
      n_totPar = edim + aux_dim*NAtoms
      ALLOCATE (n_pWrap(edim))
      n_pWrap = 0
      zdmin = 1d10
      zdmax = 0d0

      DO i = 1, NAtoms
         IF (z_priorType(i) == 0) THEN
            IF (z_prior(i, 1) < zdmin) zdmin = z_prior(i, 1)
            IF (z_prior(i, 1) > zdmax) zdmax = z_prior(i, 1)
         ELSEIF (z_priorType(i) == 1 .or. z_priorType(i) == 2 .or. z_priorType(i) == 8) THEN
            IF (z_prior(i, 1) < zdmin) zdmin = z_prior(i, 1)
            IF (z_prior(i, 2) > zdmax) zdmax = z_prior(i, 2)
         ELSEIF (z_priorType(i) == 3 .or. z_priorType(i) == 4) THEN
            IF (max(0.01d0, z_prior(i, 1) - 5d0*z_prior(i, 2)) < zdmin) zdmin = max(0.01d0, z_prior(i, 1) - 5d0*z_prior(i, 2))
            IF (z_prior(i, 1) + 5d0*z_prior(i, 2) > zdmax) zdmax = z_prior(i, 1) + 5d0*z_prior(i, 2)
         ELSE
            if (myID == 0) then
               WRITE (*, *) "ERROR: Can not use a prior other than delta, unfiorm, log uniform, Gaussian, &
               &                                        lognormal or mass function for the cluster redshift. Aborting."
            end if
            Initialise = 1
            return
         END IF

         ! uniform sampling in triangle check

         IF ((Geo_PriorType(i, 1) == 9 .and. Geo_PriorType(i, 2) /= 9) .or. &
             (Geo_PriorType(i, 1) /= 9 .and. Geo_PriorType(i, 2) == 9)) THEN
            if (myID == 0) then
               WRITE (*, *) "ERROR: Geo_PriorType should be set to 9 for both x and y positions if &
               &                                                uniform sampling in a triangle is desired"
            end if
            Initialise = 1
            return
         ELSEIF (Geo_PriorType(i, 1) == 9 .and. Geo_PriorType(i, 2) == 9) THEN
            DO j = 1, 2
               DO k = j + 1, 3
                  IF (Geo_Prior(i, j, 1) == Geo_Prior(i, k, 1) .and. Geo_Prior(i, j, 2) == Geo_Prior(i, k, 2)) THEN
                     if (myID == 0) WRITE (*, *) "ERROR: Two vertices of the triangle can not be same"
                     Initialise = 1
                     return
                  END IF
               END DO
            END DO

            DO j = 1, 2
               idx = minloc(Geo_Prior(i, j, 1:3))
               Geo_Tri_Prior(i, j, 1) = Geo_Prior(i, j, idx(1))
               idx = maxloc(Geo_Prior(i, j, 1:3))
               Geo_Tri_Prior(i, j, 2) = Geo_Prior(i, j, idx(1))
            END DO
         END IF
      END DO

      IF (zdmin == zdmax) THEN
         Dn = 1
         SCZdn = 1
      ELSE
         Dn = n
         SCZdn = n
      END IF

      !WRITE(*,*)"making Dlookup done"
      ALLOCATE (lookD(Dn, 2))
      CALL makeDlookup(zdmin, zdmax)
      !WRITE(*,*)"Dlookup done"

      !allocate memory if mass function priors used for both M & z

      Mmin = 1.d90
      Mmax = 0d0

      DO i = 1, NAtoms
         IF (z_PriorType(i) == 8 .and. (Gas_PriorType(i, 1) == 8)) THEN
            IF (Gas_Prior(i, 1, 1) < 0d0) THEN
               if (myID == 0) WRITE (*, *) "ERROR: Incorrect lower limit on M200 prior. Aborting."
               initialise = 1
               return
            END IF
            IF (Gas_Prior(i, 1, 1) < Mmin) Mmin = Gas_Prior(i, 1, 1)
            IF (Gas_Prior(i, 1, 2) > Mmax) Mmax = Gas_Prior(i, 1, 2)
         END IF
      END DO

      IF (Mmax > 0d0) THEN
         ALLOCATE (lookM(n, 2), lookZ(n, n, 2))
         CALL makeMZlookup(Mmin, Mmax, zdmin, zdmax)
         !WRITE(*,*)"MZlookup done"
      END IF

! Calculating XRAYLhood0 (-log of the componet in likelihood which does not depend on the model.)

      xraytrans(1) = -float(xraynxcentre)*xraycell
      xraytrans(2) = xraycell
      xraytrans(3) = 0.0
      xraytrans(4) = -float(xraynycentre)*xraycell
      xraytrans(5) = 0.0
      xraytrans(6) = xraycell

      XRAYLhood0 = 0.0d0
      xrayCmap(1:xrayNch, 1:xraynx, 1:xrayny) = 0.d0
      xrayCpred = 0.d0

      xrayBG(1:LENx) = (xrayBG_predmax/xrayNch)* &
                       (sexpotime)*(Aeffave)*(xraycell*xraycell*sec2min*sec2min)

      ! Calculating XRAYLhood0= sum[(c_obs)_i !]

      DO i = 1, LENx
         log_Cobs_factorial = 0d0
         log_BGobs_factorial = 0d0

         IF (xrayMask(i) .EQ. 0) THEN

            IF (xrayBG_obs(i) .EQ. 0) THEN
               log_BGobs_factorial = 0d0
            ELSE
               DO j = 1, xrayBG_obs(i)
                  log_BGobs_factorial = log_BGobs_factorial + LOG(REAL(j))
               END DO
            END IF

            IF (xrayCobs(i) .EQ. 0) THEN
               log_Cobs_factorial = 0d0
            ELSE
               DO j = 1, xrayCobs(i)
                  log_Cobs_factorial = log_Cobs_factorial + LOG(REAL(j))
               END DO
            END IF

         END IF

         XRAYLhood0 = XRAYLhood0 - log_Cobs_factorial - log_BGobs_factorial
      END DO

      !if( myID == 0 ) WRITE(*,*)'     XRAYLhood0 =' , XRAYLhood0

!-----------------------------------------
! Initialise working arrays:

      ! Dynamic radius calculation
      if (rauto .eqv. .TRUE. .and. Dn == 1) then
         ! By requiring Dn = 1 we required a fixed redshift
         ! and avoid having to deal with redshift induced variation
         ! in diameter

         ! Init rmin, rlimit for fixed redshift
         angfactor = sec2rad*lookD(1, 2) ! Physical Mpc per arcsec
         r_sky_max = dble(min(xraynx, xrayny))/2*xraycell*angfactor ! radius in Mpc
         r_sky_min = 0.001*r_sky_max
         r_los_min = r_sky_min

         ! rmax needs some additional information
         ! M200 = Gas_Prior(1, 1, 2), maximum prior value
         M200_max = Gas_Prior(1, 1, 2)

         IF (znow) THEN
            rhocritz = rhocrit
         ELSE
            rhocritz = rhocritofz(zdmin)
         END IF

         ! Setting rmax to 5x R500 which is calculated as R200/1.5
         ! Estimate R200 with NFW model
         r_los_max = ((3.d0*M200_max)/(4.d0*pi*200.d0*rhocritz))**(1.d0/3.d0)/1.5d0*5.d0
         !rmax = rlimit

         write (*, *) 'Using dynamic radius limits'
         write (*, *) 'r_sky_min: ', r_sky_min, ' r_sky_max = ', r_sky_max
         write (*, *) 'r_los_min: ', r_los_min, ' r_los_max = ', r_los_max
      end if

      DO i = 1, n
         logr(i) = log10(r_los_min) + (log10(r_los_max) - log10(r_los_min))*(i - 1)/dble(n - 1)
         r(i) = 10.d0**logr(i)
      END DO

      if (maxval(r) > r_los_max) then
         write (*, *) 'r', maxval(r), 'greater than r_los_max', r_los_max, ', adjusting r_los_max to match'
         r_los_max = maxval(r)
      end if

      !if( myID == 0 ) WRITE(*,*) '         Lhood0 =', XRAYLhood0

      XRAYLhood = 0.d0

      CALL XraycountsLhood(XRAYLhood)

      NullEv = XRAYLhood

      if (myID == 0) then
         WRITE (*, *) '  log(null-evidence) =', NullEv
         WRITE (*, *)
         WRITE (*, *) '  Calling MultiNest .........'
         WRITE (*, *)
         WRITE (*, *)
      end if

   END function Initialise

!=======================================================================

   SUBROUTINE Welcome

      WRITE (*, '(a)') '****************************************************************'
      WRITE (*, '(a)')
      WRITE (*, '(a)') '  BAYES-X v 1.3'
      WRITE (*, '(a)')
      WRITE (*, '(a)') '  Bayesian inference tool for the analysis of X-ray'
      WRITE (*, '(a)') '  observations of galaxy clusters'
      WRITE (*, '(a)')
      WRITE (*, '(a)') '  M. Olamaie et al. '
      WRITE (*, '(a)') '  September 2014 '
      WRITE (*, '(a)')
      WRITE (*, '(a)') '  Updated by R. Cox'
      WRITE (*, '(a)') '  2022-2023'
      WRITE (*, '(a)')
      WRITE (*, '(a)') '****************************************************************'
      WRITE (*, '(a)')

      RETURN
   END SUBROUTINE Welcome

!=======================================================================

END MODULE like
