MODULE massfunction

   USE cosmology
   USE constants
   USE params

   REAL*8                         ::  zglob, MM, Mglob, Minfglob

CONTAINS

!=======================================================================

   FUNCTION PSfunc(M, zloc)

      !The Evrard approximation to the halo mass function, as described in
      !Allen et al 2002. Function returns dn/dMdz in h^4 M_sun^-1 Mpc^-3. Flat
      !cosmology only at present.

      IMPLICIT NONE

      REAL*8                         ::  M, zloc, PSfunc
      REAL*8                         ::  Omz, g0, gz, sigma8z, sigmaz, rhobar, rhobarz, R
      REAL*8                         ::  BigGamma, LittleGamma, x, A, B, c, eps, f, dndM

!-----------------------------------------------------------------------

      ! First find density parameter at redshift zloc:

      Omz = Omofz(Om, zloc)

      ! Now find growth factors at redshift 0 and zloc:

      g0 = growth(0.d0)
      gz = growth(zloc)

      ! Now find sigma8 at redshift zloc:

      sigma8z = sigma8*(gz/g0)/(1.+zloc)

      ! Now need the radius which corresponds to M...

      rhobar = Om*rhocrit
      R = (3.*M/(4.*pi*rhobar))**(1./3.)

      ! ...in order to calculate shape parameters:

      BigGamma = Om*h*(2.7/2.726)*(2.7/2.726)*exp(-Ob - sqrt(h/0.5)*Ob/Om)
      LittleGamma = (0.3*BigGamma + 0.2)*(2.92 + log10(R*h/8.))

      ! So now can calculate pieces to go in Evrard formula:

      rhobarz = rhobar*(1.+zloc)**3
      !Better to have comoving number density:
      !rhobarz = rhobar
      sigmaz = sigma8z*(R*h/8.)**(-LittleGamma)

      ! Interpolations for A B and eps:

      IF (mass_function == 1 .OR. mass_function == 2) THEN
         ! Evrard approximation to Press-Schechter
         IF (mass_function == 1) THEN
            x = (1.-Omz)/0.7
            A = (1.-x)*0.27 + x*0.22
            B = (1.-x)*0.65 + x*0.73
            eps = (1.-x)*3.77 + x*3.86
            ! Jenkins mass function
         ELSEIF (mass_function == 2) THEN
            A = 0.315
            B = 0.61
            eps = 3.8
         END IF
         f = A*exp(-1.*abs(log(1./sigmaz) + B)**eps)
         !Tinker mass function
      ELSEIF (mass_function == 3) THEN
         A = 0.1d0*log10(200d0) - 0.05d0
         b = 1d0 + (log10(200d0) - 1.6d0)**(-1.5d0)
         c = 1.2d0 - (-log10(200d0) + 2.35d0)**1.6d0
         eps = 1.43d0 + (log10(200d0) - 2.3d0)**1.5d0
         f = A*(((sigmaz/b)**(-eps)) + 1d0)*exp(-c/(sigmaz**2d0))
      ELSE
         if (myID == 0) WRITE (*, *) "ERROR: Incorrect value set for variable mass_function in the include file"
         STOP
      END IF

      dndM = LittleGamma*rhobarz*f/(3.*M*M)

      ! Done!

      PSfunc = dndM

   END FUNCTION PSfunc

!=======================================================================
!Function returns linear growth factor g(zloc), from CPT92

   FUNCTION growth(zloc)

      IMPLICIT NONE

      REAL*8                         ::  growth, zloc, Omz

      Omz = Omofz(Om, zloc)
      growth = (5./2.)*Omz/(1./70.+209.*Omz/140.-Omz**2/140.+Omz**(4./7.))

   END FUNCTION growth

!=======================================================================
!Function returns dn(logM,zloc)/dlogM in h^3 Mpc^-3

   FUNCTION PSlogfunc(logM, zloc)

      IMPLICIT NONE

      REAL*8                         ::  PSlogfunc, logM, zloc
      REAL*8                         ::  M

      M = exp(logM)
      PSlogfunc = M*PSfunc(M, zloc)

   END FUNCTION PSlogfunc

!=======================================================================
!Integrate over masses to get n(>M,zloc) in h^3 Mpc^-3

   SUBROUTINE ngtMofz(M, Minfinity, zloc, Result)

      IMPLICIT NONE

      REAL*8                         ::  M, zloc, Result
      REAL*8                         ::  Minfinity, check1, check2
      !PARAMETER(Minfinity=1.d17)
      REAL*8, PARAMETER            ::  eps = 1d-3

      zglob = zloc
      check1 = ngtMofzIntegrand(log(M))
      check2 = ngtMofzIntegrand(log(Minfinity))
      IF (check1 < 1d-12 .AND. check2 < 1d-12) THEN
         Result = 0.
      ELSE
         CALL qtrap(ngtMofzIntegrand, log(M), log(Minfinity), eps, Result)
      END IF

   END SUBROUTINE ngtMofz

!-----------------------------------------------------------------------

   FUNCTION ngtMofzIntegrand(logM)

      IMPLICIT NONE

      REAL*8                         ::  ngtMofzIntegrand, logM

      ngtMofzIntegrand = PSlogfunc(logM, zglob)

   END FUNCTION ngtMofzIntegrand

!=======================================================================
!Differential  dN(M,<zloc)/dlogM in Msun^-1 sr^-1

   SUBROUTINE dNltzdM(M, zloc, Result)

      IMPLICIT NONE

      REAL*8                         ::  M, zloc, Result
      REAL*8, PARAMETER             ::  eps = 1d-3

      MM = log(M)
      CALL qtrap(dNltzdMIntegrand, 0.d0, zloc, eps, Result)

   END SUBROUTINE dNltzdM

!-----------------------------------------------------------------------

   FUNCTION dNltzdMIntegrand(zloc)

      IMPLICIT NONE

      REAL*8                         ::  zloc, dNltzdMIntegrand

      dNltzdMIntegrand = PSlogfunc(MM, zloc)
      dNltzdMIntegrand = dNltzdMIntegrand*dVdzcomoving(zloc)

   END FUNCTION dNltzdMIntegrand

!=======================================================================
!Differential  dN(>M,zloc)/dz in sr^-1

   SUBROUTINE dNgtMdz(M, Minfinity, zloc, Result)

      IMPLICIT NONE

      REAL*8                         ::  M, Minfinity, zloc, Result

      CALL ngtMofz(M, Minfinity, zloc, Result)

      Result = Result*dVdzcomoving(zloc)

   END SUBROUTINE dNgtMdz

!=======================================================================
!Integral  N(>M,<zloc) in sr^-1

   SUBROUTINE NgtMltz(M, Minfinity, zloc, Result)

      IMPLICIT NONE

      REAL*8                         ::  M, Minfinity, zloc, Result
      REAL*8, PARAMETER             ::  eps = 1d-3

      MM = M
      Minfglob = Minfinity

      CALL qtrap(NgtMltzIntegrand, 0.d0, zloc, eps, Result)

   END SUBROUTINE NgtMltz

!-----------------------------------------------------------------------

   FUNCTION NgtMltzIntegrand(zloc)

      IMPLICIT NONE

      REAL*8                         :: NgtMltzIntegrand, zloc
      REAL*8                         :: Minfinity

      Minfinity = Minfglob

      CALL dNgtMdz(MM, Minfinity, zloc, NgtMltzIntegrand)

   END FUNCTION NgtMltzIntegrand

!-----------------------------------------------------------------------
   SUBROUTINE makeMZlookup(Mmin, Mmax, zmin, zmax)

      IMPLICIT NONE

      INTEGER                         :: i, j, k
      REAL*8                          ::Mmin, Mmax, Minc, M, zmin, zmax, zinc, zz
      REAL*8                         :: dNdMdz(n, n)
      REAL*8                          ::normM, normz

!        Mmin=Mass_Prior(1,2,1) !smallest mass (prior)
!              Mmax=Mass_Prior(1,2,2) !largest mass (prior)
!              zmin=zdmin !smallest redshift (prior)
!              zmax=zdmax !largest redshift (prior)

      !set up the step size
      !n is the no. of steps, generally set to 512
      Minc = (Mmax - Mmin)/n
      zinc = (zmax - zmin)/n

      DO i = 1, n
         !look-up table for M & its cumulative probability
         !lookM(i,1)=Mmin+Minc/2.+Minc*(i-1.)
         IF (i == n) THEN
            lookM(i, 1) = Mmax
         ELSE
            lookM(i, 1) = 10.**(log10(Mmin) + (log10(Mmax) - log10(Mmin))*(i - 1.)/(n - 1.))
         END IF

         !calculate the cumulative Pr(M)
         IF (i == 1) THEN
            lookM(i, 2) = cumPrM(log(Mmin), log(lookM(i, 1)))
         ELSE
            lookM(i, 2) = lookM(i - 1, 2) + cumPrM(log(lookM(i - 1, 1)), log(lookM(i, 1)))
         END IF

         !now calculate the cumulative probability of z given M
         !set Mglob so that M is fixed
         Mglob = log(lookM(i, 1))

         !iterate over the z values from zmin to zmax
         DO j = 1, n
            !look-up table for M & its cumulative probability
            !lookZ(i,j,1)=zmin+zinc/2.+zinc*(j-1.)
            IF (j == n) THEN
               lookZ(i, j, 1) = zmax
            ELSE
               lookZ(i, j, 1) = 10.**(log10(zmin) + (log10(zmax) - log10(zmin))*(j - 1.)/(n - 1.))
            END IF

            !calculate the cumulative Pr(z)
            IF (j == 1) THEN
               IF (zmin == zmax) THEN
                  lookZ(i, j, 2) = 0d0
               ELSE
                  lookZ(i, j, 2) = cumPrz(zmin, lookZ(i, j, 1))
               END IF
            ELSE
               IF (zmin == zmax) THEN
                  lookZ(i, j, 2) = 1d0
               ELSE
                  lookZ(i, j, 2) = lookZ(i, j - 1, 2) + cumPrz(lookZ(i, j - 1, 1), lookZ(i, j, 1))
               END IF
            END IF
         END DO

         !normalizing constants for the probabilities
         normZ = lookZ(i, n, 2) + cumPrz(lookZ(i, n, 1), zmax)

         !normalize the cumulative probabilities of z given M
         lookZ(i, 1:n, 2) = lookZ(i, 1:n, 2)/normz
      END DO

      !normalizing constants for the probabilities
      normM = lookM(n, 2) + cumPrM(log(lookM(n, 1)), log(Mmax))

      !normalize the cumulative probabilities of M
      lookM(1:n, 2) = lookM(1:n, 2)/normM

   END SUBROUTINE makeMZlookup

!-----------------------------------------------------------------------
!dNdMdz returns the dN/dlogMdz
   FUNCTION dNdMdz(logM, zz)

      IMPLICIT NONE

      REAL*8                         :: dNdMdz, logM, zz

      dNdMdz = PSlogfunc(logM, zz)*dVdzcomoving(zz)

   END FUNCTION dNdMdz

!-----------------------------------------------------------------------
!dNdMdz returns the dN/dlogMdz with the M value taken as the global one & z passed
   FUNCTION dNdMdz1(zz)

      IMPLICIT NONE

      REAL*8                         :: dNdMdz1, zz, logM

      logM = Mglob
      dNdMdz1 = dNdMdz(logM, zz)

   END FUNCTION dNdMdz1

!-----------------------------------------------------------------------
!Prz returns the un-normalized probability of z given M which is equal to dN/DMdz
!z is passed as an argument & M is taken to be the global variable Mglob
   FUNCTION Prz(zz)

      IMPLICIT NONE

      REAL*8                         :: Prz, zz, logM

      logM = Mglob
      Prz = dNdMdz(logM, zz)

   END FUNCTION Prz

!-----------------------------------------------------------------------
!PrM returns the un-normalized cumulative probability of z given M
   FUNCTION cumPrz(zmin, zz)

      IMPLICIT NONE

      REAL*8                         :: zmin, zz, cumPrz
      REAL*8, PARAMETER             :: eps = 1d-3

      CALL qtrap(Prz, zmin, zz, eps, cumPrz)

   END FUNCTION cumPrz

!-----------------------------------------------------------------------
!PrM returns the un-normalized probability of M by marginalizing dN/DMdz over z
   FUNCTION PrM(logM)

      IMPLICIT NONE

      REAL*8                         :: logM, zmin, zmax, PrM
      REAL*8, PARAMETER             :: eps = 1d-3

      !get the min & max redshift values from the prior
      zmin = zdmin
      zmax = zdmax

      !set the global M value, so that qtrap calls a function with takes only 1 variable
      Mglob = logM

      IF (zmin == zmax) THEN
         PrM = dNdMdz1(zmin)
      ELSE
         CALL qtrap(dNdMdz1, zmin, zmax, eps, PrM)
      END IF

   END FUNCTION PrM

!-----------------------------------------------------------------------
!PrM returns the un-normalized cumulative probability of M
   FUNCTION cumPrM(logMmin, logM)

      IMPLICIT NONE

      REAL*8                         :: logMmin, logM, cumPrM
      REAL*8, PARAMETER             :: eps = 1d-3

      CALL qtrap(PrM, logMmin, logM, eps, cumPrM)

   END FUNCTION cumPrM

!-----------------------------------------------------------------------

   SUBROUTINE lookUpMZ(cPM, cPZ, M, zz)

      IMPLICIT NONE

      REAL*8                          :: cPM !cumulative probability of M
      REAL*8                          :: cPZ !cumulative probability of z given M
      REAL*8                          :: M !M corresponding to cPM
      REAL*8                          :: zz !z corresponding to cPZ
      INTEGER                         :: i, j, k
      REAL*8                          :: zz1, zz2

      !Lookup the Mass table
      i = binSearch(n, lookM(:, 2), cPM)

      !Lookup the z table for M>M_req
      j = binSearch(n, lookZ(i, :, 2), cPZ)

      !Lookup the z table for M<M_req
      k = binSearch(n, lookZ(i - 1, :, 2), cPZ)

      !sanity check
      IF (i == 0 .or. j == 0 .or. k == 0 .or. i > n .or. j > n .or. k > n) THEN
         if (myID == 0) WRITE (*, *) "ERROR: problem in lookUpMZ, value to look for is not within the table bounds"
         stop
      END IF

      !interpolation to get the z values from z tables
      zz1 = lookZ(i, j - 1, 1) + (lookZ(i, j, 1) - lookZ(i, j - 1, 1))*(cPZ - lookZ(i, j - 1, 2))/(lookZ(i, j, 2) - lookZ(i, j - 1, 2))
      zz2 = lookZ(i + 1, k - 1, 1) + (lookZ(i + 1, k, 1) - lookZ(i + 1, k - 1, 1))*(cPZ - lookZ(i + 1, k - 1, 2))/(lookZ(i + 1, k, 2) - lookZ(i + 1, k - 1, 2))

      !interpolation between the z values for M<M_req & M>M_req
      zz = zz1 + (zz2 - zz1)*(cPM - lookM(i - 1, 2))/(lookM(i, 2) - lookM(i - 1, 2))
      M = lookM(i - 1, 1) + (lookM(i, 1) - lookM(i - 1, 1))*(cPM - lookM(i - 1, 2))/(lookM(i, 2) - lookM(i - 1, 2))

   END SUBROUTINE lookUpMZ

!-----------------------------------------------------------------------

   SUBROUTINE lookUpM(cPM, M)

      IMPLICIT NONE

      REAL*8                          :: cPM !cumulative probability of M
      REAL*8                          :: M !M corresponding to cPM
      INTEGER                         :: i, j, k

      !Lookup the Mass table
      i = binSearch(n, lookM(:, 2), cPM)

      !sanity check
      IF (i == 0 .or. i > n) THEN
         if (myID == 0) WRITE (*, *) "ERROR: problem in lookUpM, value to look for is not within the table bounds"
         stop
      END IF

      M = lookM(i - 1, 1) + (lookM(i, 1) - lookM(i - 1, 1))*(cPM - lookM(i - 1, 2))/(lookM(i, 2) - lookM(i - 1, 2))

   END SUBROUTINE lookUpM

!-----------------------------------------------------------------------

   SUBROUTINE lookUpZ(cPZ, M, zz)

      IMPLICIT NONE

      REAL*8                          :: cPZ !cumulative probability of z given M
      REAL*8                          :: M !M corresponding to cPM
      REAL*8                          :: zz !z corresponding to cPZ
      INTEGER                         :: i, j, k
      REAL*8                          :: zz1, zz2

      !Lookup the Mass table
      i = binSearch(n, lookM(:, 1), M)

      !Lookup the z table for M>M_req
      j = binSearch(n, lookZ(i, :, 2), cPZ)

      !Lookup the z table for M<M_req
      k = binSearch(n, lookZ(i - 1, :, 2), cPZ)

      !sanity check
      IF (i == 0 .or. j == 0 .or. k == 0 .or. i > n .or. j > n .or. k > n) THEN
         if (myID == 0) WRITE (*, *) "ERROR: problem in lookUpZ, value to look for is not within the table bounds"
         stop
      END IF

      !interpolation to get the z values from z tables
      zz1 = lookZ(i, j - 1, 1) + (lookZ(i, j, 1) - lookZ(i, j - 1, 1))*(cPZ - lookZ(i, j - 1, 2))/(lookZ(i, j, 2) - lookZ(i, j - 1, 2))
      zz2 = lookZ(i + 1, k - 1, 1) + (lookZ(i + 1, k, 1) - lookZ(i + 1, k - 1, 1))*(cPZ - lookZ(i + 1, k - 1, 2))/(lookZ(i + 1, k, 2) - lookZ(i + 1, k - 1, 2))

      !interpolation between the z values for M<M_req & M>M_req
      zz = zz1 + (zz2 - zz1)*(M - lookM(i - 1, 1))/(lookM(i, 1) - lookM(i - 1, 1))

   END SUBROUTINE lookUpZ

!-----------------------------------------------------------------------

!=======================================================================

END MODULE massfunction
