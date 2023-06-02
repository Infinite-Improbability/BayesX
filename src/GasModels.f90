MODULE GasModels
   USE params
   USE constants
   USE utilities
   USE matrix_utils
   USE cosmology
   USE ReadWrite

   IMPLICIT NONE
   REAL*8                       ::  lmin, lmax, gd1, gd2
   real*8, private              :: Gamma0, GammaR, T0_poly

CONTAINS

!===========================================================
   SUBROUTINE MakeXrayGasDistributions(k, flag)

      IMPLICIT NONE

      INTEGER                       ::  i, j, k, m, flag

      REAL*8                       ::  yc, index, index_g
      REAL*8                       ::  rs, ps, cc, T0, logT0, rlimit1, rmin_fraction
      REAL*8                       ::  prefactor, prefactor2
      REAL*8                       ::  thetaE
      REAL*8, PARAMETER           ::  eps = 1d-4
      REAL*8                       ::  result
      REAL*8                       ::  angfactor

      real*8 :: ne500_poly, ne2500_poly, ne200_poly, ne_rx
!-----------------------------------------------------------------------

      ! Get redshift and a corresponding critical density
      IF (znow) THEN
         rhocritz = rhocrit
      ELSE
         rhocritz = rhocritofz(z(k))
      END IF

      ! Get angular diameter distance at current redshift
      CALL lookUp1D(Dn, z(k), lookD(:, 1), lookD(:, 2), D)

      ! Select gass model
      ! Model 1 is the NFW-GNFW model
      IF (GasModel == 1) THEN
         ! Load in current values from the priors
         MT200_DM = GasPars(k, 1)   !M_sun
         fg200_DM = GasPars(k, 2)
         a_GNFW = GasPars(k, 3)
         b_GNFW = GasPars(k, 4)
         c_GNFW = GasPars(k, 5)
         c500_GNFW = GasPars(k, 6)
         rmin_fraction = GasPars(k, 11)

         ! Sanity check on priors
         IF (fg200_DM .LT. 0.0 .OR. MT200_DM .LT. 0.0 .OR. a_GNFW .LE. 0.0 .OR. c500_GNFW .LE. 0.0 .OR. (b_GNFW - c_GNFW) .LE. 0.0 .OR. rmin_fraction .LE. 0) THEN
            flag = 1
            RETURN
         END IF

         ! Calculate c200 from mass
         ! This is based on N-body simulations
         !Neto et~al. 2007 for relaxed clusters
         c200_DM = 5.26d0*(((MT200_DM*h)/1.d14)**(-0.1d0))*(1.d0/(1.d0 + z(k)))

         ! Calculate gas mass at R200
         Mg200_DM = MT200_DM*fg200_DM     !M_sun

         ! Sanity check - null run
         IF (Mg200_DM == 0.d0) THEN
            flag = 2
            RETURN
         END IF

         ! Calculate radius from M200
         ! TODO: I think this is obtained from the NFW equation. Verify
         r200_DM = ((3.d0*MT200_DM)/(4.d0*pi*200.d0*rhocritz))**(1.d0/3.d0)   !Mpc

         ! Calculate the NFW scale radius
         rs_DM = r200_DM/c200_DM !Mpc

         ! Adjust integration limits
         ! This should, I hope, mask r_min in subsequent equations.
         r_min = rs_DM * rmin_fraction

         ! Sanity check
         if (r_min < minval(r)) then
            flag = 3
            return
         end if

         ! Calculate the dark matter density at R200
         ! TODO: Where's this equation from exactly?
         rhos_DM = (200.d0/3.d0)*((r200_DM/rs_DM)**3.d0)* &
            (rhocritz/(DLOG(1.d0 + r200_DM/rs_DM) - (1.d0/(1.d0 + rs_DM/r200_DM))))   !M_sun Mpc-3

         ! Approximate r500 from r200
         ! Also get appropriate c500 (of dark matter, do not confuse with the c500 used in GNFW)
         ! TODO: Try actually calculating it
         r500_DM = r200_DM/1.5d0                 !Mpc.
         c500_DM = r500_DM/rs_DM

         ! Set the GNFW scale radius
         rp_GNFW = r500_DM/c500_GNFW    !Mpc

         ! Calculate Pei, a normalisation coefficent for pressure profile
         ! TODO: Verify
         ! In J m-3 I guess?
         Pei_GNFW = ((mu_m/mu_e)*(G*rhos_DM*rs_DM*rs_DM*rs_DM)*Mg200_DM)/ &
            DM_GNFWgasVol(r200_DM, rs_DM, rp_GNFW, a_GNFW, b_GNFW, c_GNFW)

         ! Convert Pei to keV m-3
         Pei_GNFW_keV = Pei_GNFW*(m_sun/Mpc2m)*(J2keV)       !keVm-3

         ! Sanity check
         IF (Pei_GNFW .LE. 0d0) THEN
            if (myID == 0) then
               WRITE (*, *) c200_DM, MT200_DM, fg200_DM
               WRITE (*, *) Mg200_DM, r200_DM, rs_DM
               WRITE (*, *) rhos_DM, r500_DM, rp_GNFW
               WRITE (*, *) MT500_DM, Mg500_DM, fg500_DM
               WRITE (*, *) Tg200_DM, Tg500_DM
               WRITE (*, *) 'error Pei_GNFW= ', Pei_GNFW
            end if
            flag = 1
            return
         END IF

         ! Allocate memory for arrays
         ALLOCATE (n_H(n))
         ALLOCATE (ne_nH(n), n_e(n))
         ALLOCATE (X_emiss2D(n, xrayNbin))
         ALLOCATE (logX_emiss1D(n))
         ALLOCATE (X_S1D(n))
         ALLOCATE (X_S2D(n, xrayNbin))

         ! Coefficent of the integral that gives observed surface brightness of cluster,
         ! in photons / m^2 / s / keV / arcmin^2
         ! according to eq 8 of Olamaie2015.
         xfluxsec1 = 1.0/((4.0*pi)*((1.0 + z(k))**4))

         ! Coefficent of cooling function Lambda_c
         xfluxsec2 = (3.031d-15)

         ! This shows up in eq 22 of Olamaie 2014
         ! Appears to be a coefficent of another form of the surface brightness
         xfluxsec5 = (pi*pi)/(60.0*60.0*180.0*180.0)

         ! Initialise arrays
         DO i = 1, xrayNbin
            xrayFluxCoeff(i) = 0.
            xrayFlux(i) = 0.
            xrayBinMin(i) = 0.
            xrayBinMax(i) = 0.
         END DO

         ! Initialise energy bins
         xrayBinMin(1) = xrayEmin
         xrayBinMax(1) = xrayEmin + xrayDeltaE
         DO i = 2, xrayNbin
            xrayBinMin(i) = xrayBinMin(i - 1) + xrayDeltaE
            xrayBinMax(i) = xrayBinMax(i - 1) + xrayDeltaE
         END DO

         ! Calculates the absorption to be applied to each energy bin
         ! photar is the fractional transmission
         CALL xsphab(xrayBinMin, xrayNbin, N_H_col, photar)

         ! Calculate more parmeters
         ! TODO: Verify correctness
         Rhogas500 = (mu_e/mu_m)*(1.d0/(4.d0*Pi*G))*(Pei_GNFW/rhos_DM)*(1.0d0/(rs_DM*rs_DM*rs_DM))* &
            calcneDM(r500_DM, rs_DM, rp_GNFW, a_GNFW, b_GNFW, c_GNFW)
         Tg500_DM = (4.d0*pi*0.6*m_p*G*rhos_DM)*(rs_DM*rs_DM*rs_DM)* &
            ((DLOG(1.0 + (r500_DM/rs_DM)) - (1.0/(1.0 + (rs_DM/r500_DM))))/r500_DM)* &
            (1.0 + ((r500_DM/rp_GNFW)**(a_GNFW)))* &
            (((b_GNFW*((r500_DM/rp_GNFW)**(a_GNFW))) + c_GNFW)**(-1.0)) &
            *(m_sun*Mpc2m*Mpc2m)*(J2keV)

         ! TODO: What does this do exactly?
         CALL Xray_flux_coeff(Rhogas500, Tg500_DM, n_e500, n_H500, ne_nH500, xrayBinMin, xrayBinMax, xrayNbin, xrayDeltaE, xrayFluxCoeff)

         ! Calculate some cluster properties
         n_e500 = n_e500*1.d+6
         Ke500 = Tg500_DM/(n_e500**(2.0/3.0))
         Pe500 = n_e500*Tg500_DM
         Mg500_DM = (mu_e/mu_m)*(1.d0/G)*(Pei_GNFW/rhos_DM)*(1.0d0/(rs_DM*rs_DM*rs_DM))* &
            DM_GNFWgasVol(r500_DM, rs_DM, rp_GNFW, a_GNFW, b_GNFW, c_GNFW)
         MT500_DM = (4.d0*pi/3.d0)*(500.d0*rhocritz)*(r500_DM*r500_DM*r500_DM)
         fg500_DM = Mg500_DM/MT500_DM

         ! R2500
         r2500_DM = r200_DM/3.5d0
         c2500_DM = r2500_DM/rs_DM
         Rhogas2500 = (mu_e/mu_m)*(1.d0/(4.d0*Pi*G))*(Pei_GNFW/rhos_DM)*(1.0d0/(rs_DM*rs_DM*rs_DM))* &
            calcneDM(r2500_DM, rs_DM, rp_GNFW, a_GNFW, b_GNFW, c_GNFW)
         Tg2500_DM = (4.d0*pi*0.6*m_p*G*rhos_DM)*(rs_DM*rs_DM*rs_DM)* &
            ((DLOG(1.0 + (r2500_DM/rs_DM)) - (1.0/(1.0 + (rs_DM/r2500_DM))))/r2500_DM)* &
            (1.0 + ((r2500_DM/rp_GNFW)**(a_GNFW)))* &
            (((b_GNFW*((r2500_DM/rp_GNFW)**(a_GNFW))) + c_GNFW)**(-1.0)) &
            *(m_sun*Mpc2m*Mpc2m)*(J2keV)
         CALL Xray_flux_coeff(Rhogas2500, Tg2500_DM, n_e2500, n_H2500, ne_nH2500, xrayBinMin, xrayBinMax, xrayNbin, xrayDeltaE, xrayFluxCoeff)
         n_e2500 = n_e2500*1.d+6
         Ke2500 = Tg2500_DM/(n_e2500**(2.0/3.0))
         Pe2500 = n_e2500*Tg2500_DM
         Mg2500_DM = (mu_e/mu_m)*(1.d0/G)*(Pei_GNFW/rhos_DM)*(1.0d0/(rs_DM*rs_DM*rs_DM))* &
            DM_GNFWgasVol(r2500_DM, rs_DM, rp_GNFW, a_GNFW, b_GNFW, c_GNFW)
         MT2500_DM = (4.d0*pi/3.d0)*(2500.d0*rhocritz)*(r2500_DM*r2500_DM*r2500_DM)
         fg2500_DM = Mg2500_DM/MT2500_DM

         ! R200
         Rhogas200 = (mu_e/mu_m)*(1.d0/(4.d0*Pi*G))*(Pei_GNFW/rhos_DM)*(1.0d0/(rs_DM*rs_DM*rs_DM))* &
            calcneDM(r200_DM, rs_DM, rp_GNFW, a_GNFW, b_GNFW, c_GNFW)
         Tg200_DM = (4.d0*pi*0.6*m_p*G*rhos_DM)*(rs_DM*rs_DM*rs_DM)* &
            ((DLOG(1.0 + (r200_DM/rs_DM)) - (1.0/(1.0 + (rs_DM/r200_DM))))/r200_DM)* &
            (1.0 + ((r200_DM/rp_GNFW)**(a_GNFW)))* &
            (((b_GNFW*((r200_DM/rp_GNFW)**(a_GNFW))) + c_GNFW)**(-1.0)) &
            *(m_sun*Mpc2m*Mpc2m)*(J2keV)
         CALL Xray_flux_coeff(Rhogas200, Tg200_DM, n_e200, n_H200, ne_nH200, xrayBinMin, xrayBinMax, xrayNbin, xrayDeltaE, xrayFluxCoeff)
         n_e200 = n_e200*1.d+6
         Ke200 = Tg200_DM/(n_e200**(2.0/3.0))
         Pe200 = n_e200*Tg200_DM

         ! Allocate some more memory
         ALLOCATE (rgx(13))
         ALLOCATE (rhogasx(13), n_Hx(13))
         ALLOCATE (ne_nHx(13), n_ex(13))
         ALLOCATE (Tgx(13), Kex(13), Pex(13))
         ALLOCATE (M_DMx(13), Mg_DMx(13), fg_DMx(13))

         ! ! Starting radius, in units of R500
         ! rx_incre = 0.03

         ! ! Calculate properties at points from 0.04 R500 to 1.0 R500
         ! DO m = 1, 7
         !    rx_incre = rx_incre + 0.01
         !    rgx(m) = rx_incre*r500_DM
         !    Rhogasx(m) = (mu_e/mu_m)*(1.d0/(4.d0*Pi*G))*(Pei_GNFW/rhos_DM)*(1.0d0/(rs_DM*rs_DM*rs_DM))* & !M_sunMpc^-3
         !       calcneDM(rgx(m), rs_DM, rp_GNFW, a_GNFW, b_GNFW, c_GNFW)
         !    Tgx(m) = (4.d0*pi*0.61*m_p*G*rhos_DM)*(rs_DM*rs_DM*rs_DM)* &
         !       ((DLOG(1.0 + (rgx(m)/rs_DM)) - (1.0/(1.0 + (rs_DM/rgx(m)))))/rgx(m))* &
         !       (1.0 + ((rgx(m)/rp_GNFW)**(a_GNFW)))* &
         !       (((b_GNFW*((rgx(m)/rp_GNFW)**(a_GNFW))) + c_GNFW)**(-1.0)) &
         !       *(m_sun*Mpc2m*Mpc2m)*(J2keV)
         !    CALL Xray_flux_coeff(Rhogasx(m), Tgx(m), n_ex(m), n_Hx(m), ne_nHx(m), xrayBinMin, xrayBinMax, xrayNbin, xrayDeltaE, xrayFluxCoeff)
         !    n_ex(m) = n_ex(m)*1.d+6
         !    Kex(m) = Tgx(m)/(n_ex(m)**(2.0/3.0))
         !    Pex(m) = n_ex(m)*Tgx(m)
         !    M_DMx(m) = calcDMmass(rs_DM, rhos_DM, rgx(m))
         !    Mg_DMx(m) = (mu_e/mu_m)*(1.d0/G)*(Pei_GNFW_keV*Mpc2m/(rhos_DM*m_sun*J2keV))*(1.0d0/(rs_DM*rs_DM*rs_DM))* &
         !       DM_GNFWgasVol(rgx(m), rs_DM, rp_GNFW, a_GNFW, b_GNFW, c_GNFW)         !M_sun
         !    fg_DMx(m) = Mg_DMx(m)/M_DMx(m)
         ! END DO

         ! ! Calculate properties at points from 1.05 R500 to 1.3 R500
         ! rx_incre = 0.1
         ! DO m = 8, 13
         !    rx_incre = rx_incre + 0.05
         !    rgx(m) = rx_incre*r500_DM
         !    Rhogasx(m) = (mu_e/mu_m)*(1.d0/(4.d0*Pi*G))*(Pei_GNFW/rhos_DM)*(1.0d0/(rs_DM*rs_DM*rs_DM))* & !M_sunMpc^-3
         !       calcneDM(rgx(m), rs_DM, rp_GNFW, a_GNFW, b_GNFW, c_GNFW)
         !    Tgx(m) = (4.d0*pi*0.61*m_p*G*rhos_DM)*(rs_DM*rs_DM*rs_DM)* &
         !       ((DLOG(1.0 + (rgx(m)/rs_DM)) - (1.0/(1.0 + (rs_DM/rgx(m)))))/rgx(m))* &
         !       (1.0 + ((rgx(m)/rp_GNFW)**(a_GNFW)))* &
         !       (((b_GNFW*((rgx(m)/rp_GNFW)**(a_GNFW))) + c_GNFW)**(-1.0)) &
         !       *(m_sun*Mpc2m*Mpc2m)*(J2keV)
         !    CALL Xray_flux_coeff(Rhogasx(m), Tgx(m), n_ex(m), n_Hx(m), ne_nHx(m), xrayBinMin, xrayBinMax, xrayNbin, xrayDeltaE, xrayFluxCoeff)
         !    n_ex(m) = n_ex(m)*1.d+6
         !    Kex(m) = Tgx(m)/(n_ex(m)**(2.0/3.0))
         !    Pex(m) = n_ex(m)*Tgx(m)
         !    Mg_DMx(m) = (mu_e/mu_m)*(1.d0/G)*(Pei_GNFW_keV*Mpc2m/(rhos_DM*m_sun*J2keV))*(1.0d0/(rs_DM*rs_DM*rs_DM))* &
         !       DM_GNFWgasVol(rgx(m), rs_DM, rp_GNFW, a_GNFW, b_GNFW, c_GNFW)
         !    M_DMx(m) = calcDMmass(rs_DM, rhos_DM, rgx(m))
         !    fg_DMx(m) = Mg_DMx(m)/M_DMx(m)
         ! END DO

         write(*,*) ''
         write(*,*) '**********************************'
         write(*,*) ''

         ! Recalculate gas density, temperature and Xray_flux_coeff at points determined by logr in like.f90
         ! The Xray_flux_coeff for each point is stored in a 2D array of (radius_index, energy_bin)
         DO m = 1, n
            Rhogas(m) = (mu_e/mu_m)*(1.d0/(4.d0*Pi*G))*(Pei_GNFW/rhos_DM)*(1.0d0/(rs_DM*rs_DM*rs_DM))* & !M_sunMpc^-3
               calcneDM(r(m), rs_DM, rp_GNFW, a_GNFW, b_GNFW, c_GNFW)
            T(m) = (4.d0*pi*0.61*m_p*G*rhos_DM)*(rs_DM*rs_DM*rs_DM)* &
               ((DLOG(1.0 + (r(m)/rs_DM)) - (1.0/(1.0 + (rs_DM/r(m)))))/r(m))* &
               (1.0 + ((r(m)/rp_GNFW)**(a_GNFW)))* &
               (((b_GNFW*((r(m)/rp_GNFW)**(a_GNFW))) + c_GNFW)**(-1.0)) &
               *(m_sun*Mpc2m*Mpc2m)*(J2keV)

            CALL Xray_flux_coeff(Rhogas(m), T(m), n_e(m), n_H(m), ne_nH(m), xrayBinMin, xrayBinMax, xrayNbin, xrayDeltaE, xrayFluxCoeff)
            X_emiss2D(m, 1:xrayNbin) = xrayFluxCoeff(1:xrayNbin)

            write(*,*) m, r(m), Rhogas(m), T(m),  xrayFluxCoeff(m)
         END DO

         ! Clear up some memory
         ! TODO: Can we deallocate more?
         DEALLOCATE (n_H, n_e, ne_nH)

         ! For each energy bin E calculate S_x(s, E) (eq. 22)??
         DO i = 1, xrayNbin

            ! Looping over points in radius
            DO m = 1, n
               ! Why are we applying the cooling function coefficent now?
               ! This array will be used by XraySintegrand via Xrayemissfunc1
               logX_emiss1D(m) = phlog10(X_emiss2D(m, i)*xfluxsec2*m2cm*m2cm*m2cm*Mpc2m*Mpc2m*Mpc2m)
            END DO

            ! Looping over points in radius
            DO m = 1, n
               ! We're finding the limits on our integration radius based on our current radius.
               ! Given that the formula seems to just be Pythagorus, we must be treating the max
               ! radius as the hypotenuse, which implies uu is either the radius in the sky plane
               ! or the radius along the line of sight.
               ! As rlimit1 is used as the limits of the surface brightness equation (eq 22)
               ! which is integrated along the line of sight we can assume it must be LOS radius,
               ! so uu is the sky-plane radius.
               uu = r(m)
               rlimit1 = sqrt(max(r_integration_max*r_integration_max - uu*uu, 0.d0))

               ! Actually compute the integral for current sky-plane radius (index m)
               ! and energy (index i)
               IF (rlimit1 > 0d0) THEN
                  CALL qtrap(XraySintegrand, -rlimit1, rlimit1, eps, X_S1D(m))
               END IF

               ! Change the units and save
               ! This is predicted count rate in some region (spatial, energy, time, more dimensions?)
               ! TODO: Figure out the exact units
               X_S2D(m, i) = X_S1D(m)/(Mpc2m*Mpc2m*m2cm*m2cm)
            END DO
         END DO

         ! Apply coefficents and photoelectric absorption from foreground gases
         DO i = 1, xrayNbin
            DO m = 1, n
               X_S2D(m, i) = X_S2D(m, i)*xfluxsec1*xfluxsec5*photar(i)
            END DO
         END DO

         ! Free up some memory
         DEALLOCATE (X_emiss2D, logX_emiss1D, X_S1D)

         !=======================================================================
         ! Calculating the telescope predicted output counts for each pixel in each channel

         ALLOCATE (predX_S2D(n, xrayNch))

         ! Apply response function to update predicted count rate
         DO m = 1, n
            DO i = 1, xrayNch
               predX_S2D(m, i) = 0.0
               DO j = 1, xrayNbin
                  predX_S2D(m, i) = predX_S2D(m, i) + TRM(j, i)*X_S2D(m, j)
               END DO
            END DO
         END DO

         ! Free up memory
         DEALLOCATE (X_S2D)

         ! Load up the details of positioning
         angfactor = sec2rad*D ! Physical Mpc per arcsec
         xrayx0 = GeoPars(k, 1) ! x coordinate of cluster center in arcsec
         xrayy0 = GeoPars(k, 2) ! y coordinate of cluster center in arcsec

         ! Initalise array values
         xrayCmap = 0d0 ! Predicted counts from source
         xrayCpred = 0d0 ! Predicted counts from source and background

         ! Iterate over every pixel (cell) in the observation
         DO xrayxpix = 1, xraynx
            DO xrayypix = 1, xrayny

               ! Copied from like.f90 for clarity, this is not a disabled part of GasModels.f90
               ! xraytrans(1) = -float(xraynxcentre)*xraycell
               ! xraytrans(2) = xraycell
               ! xraytrans(3) = 0.0
               ! xraytrans(4) = -float(xraynycentre)*xraycell
               ! xraytrans(5) = 0.0
               ! xraytrans(6) = xraycell

               ! Set position
               xrayxx = xraytrans(1) + xrayxpix*xraytrans(2) ! x position in arcsecs
               xrayyy = xraytrans(4) + xrayypix*xraytrans(6) ! y position in arcsecs
               xraydx(1) = xrayxx - xrayx0 ! x position relative to center in arcsecs
               xraydx(2) = xrayyy - xrayy0 ! y position relative to center in arcsecs
               xrayr = SQRT(xraydx(1)*xraydx(1) + xraydx(2)*xraydx(2))*angfactor ! distance from center in arcsecs

               ! For each energy channel get counts by interpolating over predicted count rate
               DO i = 1, xrayNch
                  IF (xrayr < r_min) THEN
                     ! If radius is below the minimum, pretend we are at the minimum
                     call interp1d_even(predX_S2D(1:n, i), logr, n, phlog10(r_min), result)
                     xrayCmap(i, xrayxpix, xrayypix) = result
                  ELSEIF (xrayr > r_integration_max) then ! TODO: Is just > safe?
                     ! If radius is outside the sky set counts to zero
                     ! If r_sky_max < r_integration_max then this should still be impossible
                     ! But if that doesn't hold...
                     xrayCmap(i, xrayxpix, xrayypix) = 0.
                  ELSE
                     ! If we don't trigger any of the special cases just do the interpolation
                     CALL interp1d_even(predX_S2D(1:n, i), logr, n, phlog10(xrayr), result)
                     xrayCmap(i, xrayxpix, xrayypix) = result
                  END IF

                  ! Turn predicted count rate into predicted count by applying exposure time and cell size
                  xrayCmap(i, xrayxpix, xrayypix) = (xrayCmap(i, xrayxpix, xrayypix))*(sexpotime)*(xraycell*xraycell*sec2min*sec2min)
               END DO
            END DO
         END DO

         ! Free up more memory
         DEALLOCATE (predX_S2D)
         xLENi = 0

         ! Add background counts to get total predicted counts
         DO xrayxpix = 1, xraynx
            DO xrayypix = 1, xrayny
               DO i = 1, xrayNch
                  xLENi = xLENi + 1
                  xrayCpred(xLENi) = xrayCmap(i, xrayxpix, xrayypix) + xrayBG(xLENi)
               END DO
            END DO
         END DO

         !Store derived parameters
         aux(k, 1) = D              !angular diameter distance in Mpc
         aux(k, 2) = rs_DM
         aux(k, 3) = rhos_DM
         aux(k, 4) = rp_GNFW
         aux(k, 5) = Pei_GNFW_keV
         aux(k, 6) = r2500_DM
         aux(k, 7) = c2500_DM
         aux(k, 8) = Mg2500_DM
         aux(k, 9) = MT2500_DM
         aux(k, 10) = fg2500_DM
         aux(k, 11) = Tg2500_DM
         aux(k, 12) = n_e2500
         aux(k, 13) = Ke2500
         aux(k, 14) = Pe2500
         aux(k, 15) = r500_DM
         aux(k, 16) = c500_DM
         aux(k, 17) = Mg500_DM
         aux(k, 18) = MT500_DM
         aux(k, 19) = fg500_DM
         aux(k, 20) = Tg500_DM
         aux(k, 21) = n_e500
         aux(k, 22) = Ke500
         aux(k, 23) = Pe500
         aux(k, 24) = r200_DM
         aux(k, 25) = c200_DM
         aux(k, 26) = Mg200_DM
         aux(k, 27) = MT200_DM
         aux(k, 28) = fg200_DM
         aux(k, 29) = Tg200_DM
         aux(k, 30) = n_e200
         aux(k, 31) = Ke200
         aux(k, 32) = Pe200

         ! Write a bunch more parameters
         ! Useful for debug
         !i = 0
         !DO m = 33, 129, 8
         !   i = i + 1
         !   aux(k, m) = rgx(i)
         !   aux(k, m + 1) = M_DMx(i)
         !   aux(k, m + 2) = Mg_DMx(i)
         !   aux(k, m + 3) = fg_DMx(i)
         !   aux(k, m + 4) = Tgx(i)
         !   aux(k, m + 5) = n_ex(i)
         !   aux(k, m + 6) = Kex(i)
         !   aux(k, m + 7) = Pex(i)
         !END DO

         ! Free up more memory
         DEALLOCATE (rgx, Tgx, n_ex, Kex, Pex)
         DEALLOCATE (rhogasx, n_Hx, ne_nHx)
         DEALLOCATE (M_DMx, Mg_DMx, fg_DMx)

         !-------------------------------------------------------------------
      ELSEIF (GasModel == 2) THEN
         MT200_DM = GasPars(k, 1)   !M_sun
         fg200_DM = GasPars(k, 2)
         a_GNFW = GasPars(k, 3)
         b_GNFW = GasPars(k, 4)
         c_GNFW = GasPars(k, 5)
         c500_GNFW = GasPars(k, 6)
         alpha_Einasto = GasPars(k, 7)
         rmin_fraction = GasPars(k, 11)

         ! Sanity check on priors
         IF (fg200_DM .LT. 0.0 .OR. MT200_DM .LT. 0.0 .OR. a_GNFW .LE. 0.0 .OR. c500_GNFW .LE. 0.0 .OR. (b_GNFW - c_GNFW) .LE. 0.0 .OR. rmin_fraction .LE. 0) THEN
            flag = 1
            RETURN
         END IF

         c200_DM = 5.26d0*(((MT200_DM*h)/1.d14)**(-0.1d0))*(1.d0/(1.d0 + z(k)))
         Mg200_DM = MT200_DM*fg200_DM     !M_sun

         !          null run
         IF (Mg200_DM == 0.d0) THEN
            flag = 2
            RETURN
         END IF

         r200_DM = ((3.d0*MT200_DM)/(4.d0*pi*200.d0*rhocritz))**(1.d0/3.d0)   !Mpc
         r_2_DM = r200_DM/c200_DM                   !Mpc

         ! Adjust integration limits
         r_min = r_2_DM * rmin_fraction

         ! Sanity check
         if (r_min < minval(r)) then
            flag = 3
            return
         end if

         Gamma_coeff1 = 3.d0/alpha_Einasto
         Gamma_coeff2 = (2.0d0/alpha_Einasto)*((r200_DM/r_2_DM)**alpha_Einasto)

         rho_2_DM = MT200_DM/( &
            4.d0*pi*DEXP(2.d0/alpha_Einasto)* &
            r_2_DM*r_2_DM*r_2_DM* &
            ((alpha_Einasto/2.0d0)**(3.d0/alpha_Einasto))* &
            (1.d0/alpha_Einasto)* &
            INCOG(Gamma_coeff1, Gamma_coeff2))

         mass_coeff_Einasto = 4.d0*pi*DEXP(2.d0/alpha_Einasto)* &
            r_2_DM*r_2_DM*r_2_DM*rho_2_DM* &
            ((alpha_Einasto/2.0d0)**(3.d0/alpha_Einasto))* &
            (1.d0/alpha_Einasto)
         r500_DM = r200_DM/1.5
         c500_DM = r500_DM/r_2_DM
         rp_GNFW = r500_DM/c500_GNFW    !Mpc
         Pei_GNFW = (mu_m/mu_e)*G*(mass_coeff_Einasto/(4.d0*pi))*Mg200_DM/ &
            EinastoDM_GNFWgasvol(r200_DM, r_2_DM, alpha_Einasto, rp_GNFW, a_GNFW, b_GNFW, c_GNFW)
         Pei_GNFW_keV = Pei_GNFW*(m_sun/Mpc2m)*(J2keV)       !keVm-3

         IF (Pei_GNFW .LE. 0d0) THEN
            WRITE (*, *) c200_DM, MT200_DM, fg200_DM
            WRITE (*, *) Mg200_DM, r200_DM, r_2_DM
            WRITE (*, *) rho_2_DM, rp_GNFW
            WRITE (*, *) 'error Pei_GNFW= ', Pei_GNFW
            STOP
         END IF

         ALLOCATE (n_H(n))
         ALLOCATE (ne_nH(n), n_e(n))
         ALLOCATE (X_emiss2D(n, xrayNbin))
         ALLOCATE (logX_emiss1D(n))
         ALLOCATE (X_S1D(n))
         ALLOCATE (X_S2D(n, xrayNbin))

         xfluxsec1 = 1.0/((4.0*pi)*((1.0 + z(k))**4))

         xfluxsec2 = (3.031d-15)

         xfluxsec5 = (pi*pi)/(60.0*60.0*180.0*180.0)

         DO i = 1, xrayNbin
            xrayFluxCoeff(i) = 0.
            xrayFlux(i) = 0.
            xrayBinMin(i) = 0.
            xrayBinMax(i) = 0.
         END DO

         xrayBinMin(1) = xrayEmin
         xrayBinMax(1) = xrayEmin + xrayDeltaE

         DO i = 2, xrayNbin
            xrayBinMin(i) = xrayBinMin(i - 1) + xrayDeltaE
            xrayBinMax(i) = xrayBinMax(i - 1) + xrayDeltaE
         END DO

         CALL xsphab(xrayBinMin, xrayNbin, N_H_col, photar)

         Gamma_coeff2 = (2.0d0/alpha_Einasto)*((r500_DM/r_2_DM)**alpha_Einasto)
         Rhogas500 = (mu_e/mu_m)*(1.d0/G)*(Pei_GNFW/mass_coeff_Einasto)* &
            (r500_DM/INCOG(Gamma_coeff1, Gamma_coeff2))* &
            ((r500_DM/rp_GNFW)**(-1.0d0*c_GNFW))* &
            ((1.0d0 + (r500_DM/rp_GNFW)**(a_GNFW))**(-1.0*(a_GNFW + b_GNFW - c_GNFW)/a_GNFW))* &
            ((b_GNFW*((r500_DM/rp_GNFW)**(a_GNFW))) + c_GNFW)
         Tg500_DM = (0.6d0*m_p*G)*(mass_coeff_Einasto)* &
            (INCOG(Gamma_coeff1, Gamma_coeff2)/r500_DM)* &
            (1.0d0 + ((r500_DM/rp_GNFW)**(a_GNFW)))* &
            (((b_GNFW*((r500_DM/rp_GNFW)**(a_GNFW))) + c_GNFW)**(-1.0))* &
            (m_sun*Mpc2m*Mpc2m)*(J2keV)

         CALL Xray_flux_coeff(Rhogas500, Tg500_DM, n_e500, n_H500, ne_nH500, xrayBinMin, xrayBinMax, &
            xrayNbin, xrayDeltaE, xrayFluxCoeff)
         n_e500 = n_e500*1.d+6
         Ke500 = Tg500_DM/(n_e500**(2.0/3.0))
         Pe500 = n_e500*Tg500_DM
         Mg500_DM = (4.d0*pi)*(mu_e/mu_m)*(1.d0/G)*(Pei_GNFW/mass_coeff_Einasto)* &
            EinastoDM_GNFWgasvol( &
            r500_DM, r_2_DM, alpha_Einasto, rp_GNFW, a_GNFW, b_GNFW, c_GNFW)
         MT500_DM = (4.d0*pi/3.d0)*(500.d0*rhocritz)*(r500_DM*r500_DM*r500_DM)
         fg500_DM = Mg500_DM/MT500_DM

         r2500_DM = r200_DM/3.5d0
         c2500_DM = r2500_DM/r_2_DM
         Gamma_coeff2 = (2.0d0/alpha_Einasto)*((r2500_DM/r_2_DM)**alpha_Einasto)
         Rhogas2500 = (mu_e/mu_m)*(1.d0/G)*(Pei_GNFW/mass_coeff_Einasto)* &
            (r2500_DM/INCOG(Gamma_coeff1, Gamma_coeff2))* &
            ((r2500_DM/rp_GNFW)**(-1.0d0*c_GNFW))* &
            ((1.0d0 + (r2500_DM/rp_GNFW)**(a_GNFW))**(-1.0*(a_GNFW + b_GNFW - c_GNFW)/a_GNFW))* &
            ((b_GNFW*((r2500_DM/rp_GNFW)**(a_GNFW))) + c_GNFW)
         Tg2500_DM = (0.6d0*m_p*G)*(mass_coeff_Einasto)* &
            (INCOG(Gamma_coeff1, Gamma_coeff2)/r2500_DM)* &
            (1.0d0 + ((r2500_DM/rp_GNFW)**(a_GNFW)))* &
            (((b_GNFW*((r2500_DM/rp_GNFW)**(a_GNFW))) + c_GNFW)**(-1.0))* &
            (m_sun*Mpc2m*Mpc2m)*(J2keV)
         CALL Xray_flux_coeff(Rhogas2500, Tg2500_DM, n_e2500, n_H2500, ne_nH2500, xrayBinMin, xrayBinMax, xrayNbin, xrayDeltaE, xrayFluxCoeff)

         n_e2500 = n_e2500*1.d+6
         Ke2500 = Tg2500_DM/(n_e2500**(2.0/3.0))
         Pe2500 = n_e2500*Tg2500_DM
         Mg2500_DM = (4.d0*pi)*(mu_e/mu_m)*(1.d0/G)*(Pei_GNFW/mass_coeff_Einasto)* &
            EinastoDM_GNFWgasvol( &
            r2500_DM, r_2_DM, alpha_Einasto, rp_GNFW, a_GNFW, b_GNFW, c_GNFW)
         MT2500_DM = (4.d0*pi/3.d0)*(2500.d0*rhocritz)*(r2500_DM*r2500_DM*r2500_DM)
         fg2500_DM = Mg2500_DM/MT2500_DM

         Gamma_coeff2 = (2.0d0/alpha_Einasto)*((r200_DM/r_2_DM)**alpha_Einasto)
         Rhogas200 = (mu_e/mu_m)*(1.d0/G)*(Pei_GNFW/mass_coeff_Einasto)* &
            (r200_DM/INCOG(Gamma_coeff1, Gamma_coeff2))* &
            ((r200_DM/rp_GNFW)**(-1.0d0*c_GNFW))* &
            ((1.0d0 + (r200_DM/rp_GNFW)**(a_GNFW))**(-1.0*(a_GNFW + b_GNFW - c_GNFW)/a_GNFW))* &
            ((b_GNFW*((r200_DM/rp_GNFW)**(a_GNFW))) + c_GNFW)
         Tg200_DM = (0.6d0*m_p*G)*(mass_coeff_Einasto)* &
            (INCOG(Gamma_coeff1, Gamma_coeff2)/r200_DM)* &
            (1.0d0 + ((r200_DM/rp_GNFW)**(a_GNFW)))* &
            (((b_GNFW*((r200_DM/rp_GNFW)**(a_GNFW))) + c_GNFW)**(-1.0))* &
            (m_sun*Mpc2m*Mpc2m)*(J2keV)
         CALL Xray_flux_coeff(Rhogas200, Tg200_DM, n_e200, n_H200, ne_nH200, xrayBinMin, xrayBinMax, xrayNbin, xrayDeltaE, xrayFluxCoeff)
         n_e200 = n_e200*1.d+6
         Ke200 = Tg200_DM/(n_e200**(2.0/3.0))
         Pe200 = n_e200*Tg200_DM

         ALLOCATE (rgx(13))
         ALLOCATE (rhogasx(13), n_Hx(13))
         ALLOCATE (ne_nHx(13), n_ex(13))
         ALLOCATE (Tgx(13), Kex(13), Pex(13))
         ALLOCATE (M_DMx(13), Mg_DMx(13), fg_DMx(13))

         rx_incre = 0.03
         DO m = 1, 7
            rx_incre = rx_incre + 0.01
            rgx(m) = rx_incre*r500_DM
            Gamma_coeff2 = (2.0d0/alpha_Einasto)*((rgx(m)/r_2_DM)**alpha_Einasto)
            Rhogasx(m) = (mu_e/mu_m)*(1.d0/G)*(Pei_GNFW/mass_coeff_Einasto)* &
               (rgx(m)/INCOG(Gamma_coeff1, Gamma_coeff2))* &
               ((rgx(m)/rp_GNFW)**(-1.0d0*c_GNFW))* &
               ((1.0d0 + (rgx(m)/rp_GNFW)**(a_GNFW))**(-1.0*(a_GNFW + b_GNFW - c_GNFW)/a_GNFW))* &
               ((b_GNFW*((rgx(m)/rp_GNFW)**(a_GNFW))) + c_GNFW)
            Tgx(m) = (0.6d0*m_p*G)*(mass_coeff_Einasto)* &
               (INCOG(Gamma_coeff1, Gamma_coeff2)/rgx(m))* &
               (1.0d0 + ((rgx(m)/rp_GNFW)**(a_GNFW)))* &
               (((b_GNFW*((rgx(m)/rp_GNFW)**(a_GNFW))) + c_GNFW)**(-1.0))* &
               (m_sun*Mpc2m*Mpc2m)*(J2keV)
            CALL Xray_flux_coeff(Rhogasx(m), Tgx(m), n_ex(m), n_Hx(m), &
               ne_nHx(m), xrayBinMin, xrayBinMax, xrayNbin, &
               xrayDeltaE, xrayFluxCoeff)
            n_ex(m) = n_ex(m)*1.d+6
            Kex(m) = Tgx(m)/(n_ex(m)**(2.0/3.0))
            Pex(m) = n_ex(m)*Tgx(m)
            Mg_DMx(m) = (4.d0*pi)*(mu_e/mu_m)*(1.d0/G)*(Pei_GNFW/mass_coeff_Einasto)* &
               EinastoDM_GNFWgasvol( &
               rgx(m), r_2_DM, alpha_Einasto, rp_GNFW, a_GNFW, b_GNFW, c_GNFW)
            M_DMx(m) = mass_coeff_Einasto*INCOG(Gamma_coeff1, Gamma_coeff2)
            fg_DMx(m) = Mg_DMx(m)/M_DMx(m)
         END DO

         rx_incre = 0.1
         DO m = 8, 13
            rx_incre = rx_incre + 0.05
            rgx(m) = rx_incre*r500_DM
            Gamma_coeff2 = (2.0d0/alpha_Einasto)*((rgx(m)/r_2_DM)**alpha_Einasto)
            Rhogasx(m) = (mu_e/mu_m)*(1.d0/G)*(Pei_GNFW/mass_coeff_Einasto)* &
               (rgx(m)/INCOG(Gamma_coeff1, Gamma_coeff2))* &
               ((rgx(m)/rp_GNFW)**(-1.0d0*c_GNFW))* &
               ((1.0d0 + (rgx(m)/rp_GNFW)**(a_GNFW))**(-1.0*(a_GNFW + b_GNFW - c_GNFW)/a_GNFW))* &
               ((b_GNFW*((rgx(m)/rp_GNFW)**(a_GNFW))) + c_GNFW)
            Tgx(m) = (0.6d0*m_p*G)*(mass_coeff_Einasto)* &
               (INCOG(Gamma_coeff1, Gamma_coeff2)/rgx(m))* &
               (1.0d0 + ((rgx(m)/rp_GNFW)**(a_GNFW)))* &
               (((b_GNFW*((rgx(m)/rp_GNFW)**(a_GNFW))) + c_GNFW)**(-1.0))* &
               (m_sun*Mpc2m*Mpc2m)*(J2keV)
            CALL Xray_flux_coeff(Rhogasx(m), Tgx(m), n_ex(m), n_Hx(m), &
               ne_nHx(m), xrayBinMin, xrayBinMax, xrayNbin, &
               xrayDeltaE, xrayFluxCoeff)
            n_ex(m) = n_ex(m)*1.d+6
            Kex(m) = Tgx(m)/(n_ex(m)**(2.0/3.0))
            Pex(m) = n_ex(m)*Tgx(m)
            Mg_DMx(m) = (4.d0*pi)*(mu_e/mu_m)*(1.d0/G)*(Pei_GNFW/mass_coeff_Einasto)* &
               EinastoDM_GNFWgasvol( &
               rgx(m), r_2_DM, alpha_Einasto, rp_GNFW, a_GNFW, b_GNFW, c_GNFW)
            M_DMx(m) = mass_coeff_Einasto*INCOG(Gamma_coeff1, Gamma_coeff2)
            fg_DMx(m) = Mg_DMx(m)/M_DMx(m)
         END DO

         DO m = 1, n
            Gamma_coeff2 = (2.0d0/alpha_Einasto)*((r(m)/r_2_DM)**alpha_Einasto)
            Rhogas(m) = (mu_e/mu_m)*(1.d0/G)*(Pei_GNFW/mass_coeff_Einasto)* &
               (r(m)/INCOG(Gamma_coeff1, Gamma_coeff2))* &
               ((r(m)/rp_GNFW)**(-1.0d0*c_GNFW))* &
               ((1.0d0 + (r(m)/rp_GNFW)**(a_GNFW))**(-1.0*(a_GNFW + b_GNFW - c_GNFW)/a_GNFW))* &
               ((b_GNFW*((r(m)/rp_GNFW)**(a_GNFW))) + c_GNFW)
            T(m) = (0.6d0*m_p*G)*(mass_coeff_Einasto)* &
               (INCOG(Gamma_coeff1, Gamma_coeff2)/r(m))* &
               (1.0d0 + ((r(m)/rp_GNFW)**(a_GNFW)))* &
               (((b_GNFW*((r(m)/rp_GNFW)**(a_GNFW))) + c_GNFW)**(-1.0))* &
               (m_sun*Mpc2m*Mpc2m)*(J2keV)
            CALL Xray_flux_coeff(Rhogas(m), T(m), n_e(m), n_H(m), ne_nH(m), &
               xrayBinMin, xrayBinMax, xrayNbin, xrayDeltaE, xrayFluxCoeff)
            X_emiss2D(m, 1:xrayNbin) = xrayFluxCoeff(1:xrayNbin)
         END DO

         DEALLOCATE (n_H, n_e, ne_nH)
         DO i = 1, xrayNbin
            DO m = 1, n
               logX_emiss1D(m) = phlog10(X_emiss2D(m, i)* &
                  xfluxsec2*m2cm*m2cm*m2cm*Mpc2m*Mpc2m*Mpc2m)
               !      write(*,*)logX_emiss1D(m)
            END DO

            DO m = 1, n
               uu = r(m)
               rlimit1 = sqrt(max(r_integration_max*r_integration_max - uu*uu, 0.d0))
               !      write(*,*)rlimit1
               IF (rlimit1 > 0d0) THEN
                  CALL qtrap(XraySintegrand, -rlimit1, rlimit1, eps, X_S1D(m))
               END IF
               X_S2D(m, i) = X_S1D(m)/(Mpc2m*Mpc2m*m2cm*m2cm)
            END DO
         END DO
         !     write(*,*)'test'
         DO i = 1, xrayNbin
            DO m = 1, n
               X_S2D(m, i) = X_S2D(m, i)*xfluxsec1*xfluxsec5*photar(i)
            END DO
         END DO
         DEALLOCATE (X_emiss2D, logX_emiss1D, X_S1D)
         !=======================================================================
         ! Calculating the telescope predicted output counts for each pixel in each channel

         ALLOCATE (predX_S2D(n, xrayNch))

         DO m = 1, n
            DO i = 1, xrayNch
               predX_S2D(m, i) = 0.0
               DO j = 1, xrayNbin
                  predX_S2D(m, i) = predX_S2D(m, i) + TRM(j, i)*X_S2D(m, j)
               END DO
            END DO
         END DO
         !  write(*,*) predX_S2D(1 ,1)
         DEALLOCATE (X_S2D)

         angfactor = sec2rad*D
         xrayx0 = GeoPars(k, 1)
         xrayy0 = GeoPars(k, 2)
         xrayCmap = 0d0
         xrayCpred = 0d0

         DO xrayxpix = 1, xraynx
            DO xrayypix = 1, xrayny
               xrayxx = xraytrans(1) + xrayxpix*xraytrans(2)
               xrayyy = xraytrans(4) + xrayypix*xraytrans(6)
               xraydx(1) = xrayxx - xrayx0
               xraydx(2) = xrayyy - xrayy0
               xrayr = SQRT(xraydx(1)*xraydx(1) + xraydx(2)*xraydx(2))*angfactor
               DO i = 1, xrayNch
                  IF (xrayr < r_min) THEN
                     CALL interp1d_even(predX_S2D(1:n, i), logr, n, phlog10(r_min), result)
                     xrayCmap(i, xrayxpix, xrayypix) = result
                  ELSEIF (xrayr > r_integration_max) then
                     xrayCmap(i, xrayxpix, xrayypix) = 0.
                  ELSE
                     !CALL interp1d(predX_S2D(1:n, i), r, n, xrayr, result)
                     CALL interp1d_even(predX_S2D(1:n, i), logr, n, phlog10(xrayr), result)
                     xrayCmap(i, xrayxpix, xrayypix) = result
                  END IF
                  xrayCmap(i, xrayxpix, xrayypix) = &
                     (xrayCmap(i, xrayxpix, xrayypix))*(sexpotime)* &
                     (xraycell*xraycell*sec2min*sec2min)
               END DO
            END DO
         END DO
         DEALLOCATE (predX_S2D)

         xLENi = 0

         DO xrayxpix = 1, xraynx
            DO xrayypix = 1, xrayny
               DO i = 1, xrayNch
                  xLENi = xLENi + 1
                  xrayCpred(xLENi) = xrayCmap(i, xrayxpix, xrayypix) + xrayBG(xLENi)
               END DO
            END DO
         END DO
         !Store derived parameters
         aux(k, 1) = D              !angular diameter distance in Mpc
         aux(k, 2) = r_2_DM
         aux(k, 3) = rho_2_DM
         aux(k, 4) = rp_GNFW
         aux(k, 5) = Pei_GNFW_keV
         aux(k, 6) = r2500_DM
         aux(k, 7) = c2500_DM
         aux(k, 8) = Mg2500_DM
         aux(k, 9) = MT2500_DM
         aux(k, 10) = fg2500_DM
         aux(k, 11) = Tg2500_DM
         aux(k, 12) = n_e2500
         aux(k, 13) = Ke2500
         aux(k, 14) = Pe2500
         aux(k, 15) = r500_DM
         aux(k, 16) = c500_DM
         aux(k, 17) = Mg500_DM
         aux(k, 18) = MT500_DM
         aux(k, 19) = fg500_DM
         aux(k, 20) = Tg500_DM
         aux(k, 21) = n_e500
         aux(k, 22) = Ke500
         aux(k, 23) = Pe500
         aux(k, 24) = r200_DM
         aux(k, 25) = c200_DM
         aux(k, 26) = Mg200_DM
         aux(k, 27) = MT200_DM
         aux(k, 28) = fg200_DM
         aux(k, 29) = Tg200_DM
         aux(k, 30) = n_e200
         aux(k, 31) = Ke200
         aux(k, 32) = Pe200

         !i = 0
         !DO m = 33, 129, 8
         !   i = i + 1
         !   aux(k, m) = rgx(i)
         !   aux(k, m + 1) = M_DMx(i)
         !   aux(k, m + 2) = Mg_DMx(i)
         !   aux(k, m + 3) = fg_DMx(i)
         !   aux(k, m + 4) = Tgx(i)
         !   aux(k, m + 5) = n_ex(i)
         !   aux(k, m + 6) = Kex(i)
         !   aux(k, m + 7) = Pex(i)
         !END DO
         DEALLOCATE (rgx, Tgx, n_ex, Kex, Pex)
         DEALLOCATE (rhogasx, n_Hx, ne_nHx)
         DEALLOCATE (M_DMx, Mg_DMx, fg_DMx)

         ! Polytropic model based on Section 4.1 of Ghirardini2019 [A&A 627, A19 (2019)]
         ! https://doi.org/10.1051/0004-6361/201834875
      ELSEIF (GasModel == 3) THEN
         MT200_DM = GasPars(k, 1)   !M_sun
         Gamma0 = GasPars(k, 8)  !dimensionless
         GammaR = GasPars(k, 9)  !dimensionless
         T0_poly = GasPars(k, 10)   !dimensionless
         rmin_fraction = GasPars(k, 11)

         ! TODO: Unit comments

         ! Sanity check on priors
         IF (fg200_DM .LT. 0.0 .OR. MT200_DM .LT. 0.0 .OR. a_GNFW .LE. 0.0 .OR. c500_GNFW .LE. 0.0 .OR. (b_GNFW - c_GNFW) .LE. 0.0 .OR. rmin_fraction .LE. 0) THEN
            flag = 1
            RETURN
         END IF

         !Neto et~al. 2007 for relaxed clusters
         c200_DM = 5.26d0*(((MT200_DM*h)/1.d14)**(-0.1d0))*(1.d0/(1.d0 + z(k)))

         !null run
         IF (MT200_DM == 0.d0) THEN
            flag = 2
            RETURN
         END IF

         ! Also used by prior models
         r200_DM = ((3.d0*MT200_DM)/(4.d0*pi*200.d0*rhocritz))**(1.d0/3.d0)   !Mpc

         ! Also used by prior models
         rs_DM = r200_DM/c200_DM ! Mpc
         r500_DM = r200_DM/1.5d0 ! Mpc
         c500_DM = r500_DM/rs_DM

         ! Adjust integration limits
         r_min = rs_DM * rmin_fraction

         ! Sanity check
         if (r_min < minval(r)) then
            flag = 3
            return
         end if

         ! Also used by prior models
         rhos_DM = (200.d0/3.d0)*((r200_DM/rs_DM)**3.d0)* &
            (rhocritz/(DLOG(1.d0 + r200_DM/rs_DM) - (1.d0/(1.d0 + rs_DM/r200_DM))))   !M_sunMpc-3

         ALLOCATE (n_H(n))
         ALLOCATE (ne_nH(n), n_e(n))
         ALLOCATE (X_emiss2D(n, xrayNbin))
         ALLOCATE (logX_emiss1D(n))
         ALLOCATE (X_S1D(n))
         ALLOCATE (X_S2D(n, xrayNbin))

         ! All used by prior models
         xfluxsec1 = 1.0/((4.0*pi)*((1.0 + z(k))**4)) ! Coefficent of integral that gives observed surface brightness of cluster, in photons / m^2 / s / keV / arcmin^2 according to eq 8 of Olamaie2015.
         xfluxsec2 = (3.031d-15)
         xfluxsec5 = (pi*pi)/(60.0*60.0*180.0*180.0)

         ! Init arrays
         DO i = 1, xrayNbin
            xrayFluxCoeff(i) = 0.
            xrayFlux(i) = 0.
            xrayBinMin(i) = 0.
            xrayBinMax(i) = 0.
         END DO

         xrayBinMin(1) = xrayEmin
         xrayBinMax(1) = xrayEmin + xrayDeltaE

         DO i = 2, xrayNbin
            xrayBinMin(i) = xrayBinMin(i - 1) + xrayDeltaE
            xrayBinMax(i) = xrayBinMax(i - 1) + xrayDeltaE
         END DO

         CALL xsphab(xrayBinMin, xrayNbin, N_H_col, photar)

         MT500_DM = (4.d0*pi/3.d0)*(500.d0*rhocritz)*(r500_DM*r500_DM*r500_DM)

         ne500_poly = polyEstimateNumberDensity(r500_DM)
         Rhogas500 = ne500_poly*mu_e/m_sun*(Mpc2m*Mpc2m*Mpc2m) ! solar masses per Mpc
         ! T500 as in eq (10) of Ghirardini et al. (2019) to match polytropic paper
         ! The last constant is actually mu/0.6, where mu is the mean molecular weight per particle
         ! I have accordingly used the molecular mass rather than SI units.
         Tg500_DM = 8.85*(MT500_DM*h/(10.**15))**(2./3.)*(rhocritz/rhocrit)**(1./3.)*(1.0078250319) ! keV

         ! TODO: n_e500 values appear to get calculated by the following, which feels odd given we've calculated it already
         CALL Xray_flux_coeff(Rhogas500, Tg500_DM, n_e500, n_H500, ne_nH500, xrayBinMin, xrayBinMax, &
            xrayNbin, xrayDeltaE, xrayFluxCoeff)
         n_e500 = n_e500*1.d+6
         Ke500 = Tg500_DM/(n_e500**(2.0/3.0))
         Pe500 = n_e500*Tg500_DM

         ! Mg500_DM = (4.d0*pi)*(mu_e/mu_m)*(1.d0/G)*(Pei_GNFW/mass_coeff_Einasto)* &
         !            EinastoDM_GNFWgasvol( &
         !            r500_DM, r_2_DM, alpha_Einasto, rp_GNFW, a_GNFW, b_GNFW, c_GNFW)
         ! fg500_DM = Mg500_DM/MT500_DM

         r2500_DM = r200_DM/3.5d0
         c2500_DM = r2500_DM/rs_DM

         ne2500_poly = polyEstimateNumberDensity(r2500_DM)
         Rhogas2500 = ne2500_poly*mu_e/m_sun*(Mpc2m*Mpc2m*Mpc2m)
         Tg2500_DM = polyTemperature(r2500_DM, ne2500_poly)

         CALL Xray_flux_coeff(Rhogas2500, Tg2500_DM, n_e2500, n_H2500, ne_nH2500, xrayBinMin, xrayBinMax, xrayNbin, xrayDeltaE, xrayFluxCoeff)
         n_e2500 = n_e2500*1.d+6

         MT2500_DM = (4.d0*pi/3.d0)*(2500.d0*rhocritz)*(r2500_DM*r2500_DM*r2500_DM)

         Ke2500 = Tg2500_DM/(n_e2500**(2.0/3.0)) ! TODO: Verify
         Pe2500 = n_e2500*Tg2500_DM

         ! Mg2500_DM = (4.d0*pi)*(mu_e/mu_m)*(1.d0/G)*(Pei_GNFW/mass_coeff_Einasto)* &
         ! EinastoDM_GNFWgasvol( &
         ! r2500_DM, r_2_DM, alpha_Einasto, rp_GNFW, a_GNFW, b_GNFW, c_GNFW)
         !fg2500_DM = Mg2500_DM/MT2500_DM

         ne200_poly = polyEstimateNumberDensity(r200_DM)
         Rhogas200 = ne200_poly*mu_e/m_sun*(Mpc2m*Mpc2m*Mpc2m)
         Tg200_DM = polyTemperature(r200_DM, ne200_poly)

         CALL Xray_flux_coeff(Rhogas200, Tg200_DM, n_e200, n_H200, ne_nH200, xrayBinMin, xrayBinMax, xrayNbin, xrayDeltaE, xrayFluxCoeff)
         n_e200 = n_e200*1.d+6

         Ke200 = Tg200_DM/(n_e200**(2.0/3.0))
         Pe200 = n_e200*Tg200_DM

         ALLOCATE (rgx(13))
         ALLOCATE (rhogasx(13), n_Hx(13))
         ALLOCATE (ne_nHx(13), n_ex(13))
         ALLOCATE (Tgx(13), Kex(13), Pex(13))
         ALLOCATE (M_DMx(13), Mg_DMx(13), fg_DMx(13))

         rx_incre = 0.03
         DO m = 1, 7
            rx_incre = rx_incre + 0.01

            rgx(m) = rx_incre*r500_DM
            ne_rx = polyEstimateNumberDensity(rgx(m))
            Rhogasx(m) = ne_rx*mu_e/m_sun*(Mpc2m*Mpc2m*Mpc2m)
            Tgx(m) = polyTemperature(rgx(m), ne_rx)

            CALL Xray_flux_coeff(Rhogasx(m), Tgx(m), n_ex(m), n_Hx(m), &
               ne_nHx(m), xrayBinMin, xrayBinMax, xrayNbin, &
               xrayDeltaE, xrayFluxCoeff)
            n_ex(m) = n_ex(m)*1.d+6
            Kex(m) = Tgx(m)/(n_ex(m)**(2.0/3.0))
            Pex(m) = n_ex(m)*Tgx(m)
            ! Mg_DMx(m) = (4.d0*pi)*(mu_e/mu_m)*(1.d0/G)*(Pei_GNFW/mass_coeff_Einasto)* &
            !             EinastoDM_GNFWgasvol( &
            !             rgx(m), r_2_DM, alpha_Einasto, rp_GNFW, a_GNFW, b_GNFW, c_GNFW)
            M_DMx(m) = calcDMmass(rs_DM, rhos_DM, rgx(m))
            ! fg_DMx(m) = Mg_DMx(m)/M_DMx(m)
         END DO

         rx_incre = 0.1
         DO m = 8, 13
            rx_incre = rx_incre + 0.05

            rgx(m) = rx_incre*r500_DM
            ne_rx = polyEstimateNumberDensity(rgx(m))
            Rhogasx(m) = ne_rx*mu_e/m_sun*(Mpc2m*Mpc2m*Mpc2m)
            Tgx(m) = polyTemperature(rgx(m), ne_rx)

            CALL Xray_flux_coeff(Rhogasx(m), Tgx(m), n_ex(m), n_Hx(m), &
               ne_nHx(m), xrayBinMin, xrayBinMax, xrayNbin, &
               xrayDeltaE, xrayFluxCoeff)
            n_ex(m) = n_ex(m)*1.d+6
            Kex(m) = Tgx(m)/(n_ex(m)**(2.0/3.0))
            Pex(m) = n_ex(m)*Tgx(m)
            ! Mg_DMx(m) = (4.d0*pi)*(mu_e/mu_m)*(1.d0/G)*(Pei_GNFW/mass_coeff_Einasto)* &
            !             EinastoDM_GNFWgasvol( &
            !             rgx(m), r_2_DM, alpha_Einasto, rp_GNFW, a_GNFW, b_GNFW, c_GNFW)
            M_DMx(m) = calcDMmass(rs_DM, rhos_DM, rgx(m))
            ! fg_DMx(m) = Mg_DMx(m)/M_DMx(m)
         END DO

         DO m = 1, n
            ne_rx = polyEstimateNumberDensity(r(m))
            Rhogas(m) = ne_rx*mu_e/m_sun*(Mpc2m*Mpc2m*Mpc2m)
            T(m) = polyTemperature(r(m), ne_rx)

            CALL Xray_flux_coeff(Rhogas(m), T(m), n_e(m), n_H(m), ne_nH(m), &
               xrayBinMin, xrayBinMax, xrayNbin, xrayDeltaE, xrayFluxCoeff)
            X_emiss2D(m, 1:xrayNbin) = xrayFluxCoeff(1:xrayNbin)
         END DO

         DEALLOCATE (n_H, n_e, ne_nH)
         DO i = 1, xrayNbin
            DO m = 1, n
               logX_emiss1D(m) = phlog10(X_emiss2D(m, i)* &
                  xfluxsec2*m2cm*m2cm*m2cm*Mpc2m*Mpc2m*Mpc2m) ! per cm^3?
               !      write(*,*)logX_emiss1D(m)
            END DO

            DO m = 1, n
               uu = r(m)
               rlimit1 = sqrt(max(r_integration_max*r_integration_max - uu*uu, 0.d0))
               !      write(*,*)rlimit1
               IF (rlimit1 > 0d0) THEN
                  CALL qtrap(XraySintegrand, -rlimit1, rlimit1, eps, X_S1D(m))
               END IF
               X_S2D(m, i) = X_S1D(m)/(Mpc2m*Mpc2m*m2cm*m2cm) ! per Mpc^2?
            END DO
         END DO
         !     write(*,*)'test'
         DO i = 1, xrayNbin
            DO m = 1, n
               X_S2D(m, i) = X_S2D(m, i)*xfluxsec1*xfluxsec5*photar(i)
            END DO
         END DO
         DEALLOCATE (X_emiss2D, logX_emiss1D, X_S1D)
         !=======================================================================
         ! Calculating the telescope predicted output counts for each pixel in each channel

         ALLOCATE (predX_S2D(n, xrayNch))

         DO m = 1, n
            DO i = 1, xrayNch
               predX_S2D(m, i) = 0.0
               DO j = 1, xrayNbin
                  predX_S2D(m, i) = predX_S2D(m, i) + TRM(j, i)*X_S2D(m, j)
               END DO
            END DO
         END DO
         !  write(*,*) predX_S2D(1 ,1)
         DEALLOCATE (X_S2D)

         angfactor = sec2rad*D
         xrayx0 = GeoPars(k, 1)
         xrayy0 = GeoPars(k, 2)
         xrayCmap = 0d0
         xrayCpred = 0d0

         DO xrayxpix = 1, xraynx
            DO xrayypix = 1, xrayny
               xrayxx = xraytrans(1) + xrayxpix*xraytrans(2)
               xrayyy = xraytrans(4) + xrayypix*xraytrans(6)
               xraydx(1) = xrayxx - xrayx0
               xraydx(2) = xrayyy - xrayy0
               xrayr = SQRT(xraydx(1)*xraydx(1) + xraydx(2)*xraydx(2))*angfactor
               DO i = 1, xrayNch
                  IF (xrayr < r_min) THEN
                     !CALL interp1d(predX_S2D(1:n, i), r, n, rmin, result)
                     CALL interp1d_even(predX_S2D(1:n, i), logr, n, phlog10(r_min), result)
                     xrayCmap(i, xrayxpix, xrayypix) = result
                  ELSEIF (xrayr > r_integration_max) THEN
                     xrayCmap(i, xrayxpix, xrayypix) = 0.
                  ELSE
                     CALL interp1d_even(predX_S2D(1:n, i), logr, n, phlog10(xrayr), result)
                     xrayCmap(i, xrayxpix, xrayypix) = result
                  END IF
                  xrayCmap(i, xrayxpix, xrayypix) = &
                     (xrayCmap(i, xrayxpix, xrayypix))*(sexpotime)* &
                     (xraycell*xraycell*sec2min*sec2min) ! multiplying by cell area in arcmin^2
               END DO
            END DO
         END DO
         DEALLOCATE (predX_S2D)

         xLENi = 0

         DO xrayxpix = 1, xraynx
            DO xrayypix = 1, xrayny
               DO i = 1, xrayNch
                  xLENi = xLENi + 1
                  xrayCpred(xLENi) = xrayCmap(i, xrayxpix, xrayypix) + xrayBG(xLENi)
               END DO
            END DO
         END DO
         !Store derived parameters
         aux(k, 1) = D              !angular diameter distance in Mpc
         aux(k, 2) = rs_DM
         aux(k, 3) = rhos_DM
         aux(k, 4) = 0!rp_GNFW
         aux(k, 5) = 0!Pei_GNFW_keV
         aux(k, 6) = r2500_DM
         aux(k, 7) = c2500_DM
         aux(k, 8) = 0!Mg2500_DM
         aux(k, 9) = MT2500_DM
         aux(k, 10) = 0!fg2500_DM
         aux(k, 11) = Tg2500_DM
         aux(k, 12) = n_e2500
         aux(k, 13) = Ke2500
         aux(k, 14) = Pe2500
         aux(k, 15) = r500_DM
         aux(k, 16) = c500_DM
         aux(k, 17) = 0!Mg500_DM
         aux(k, 18) = MT500_DM
         aux(k, 19) = 0!fg500_DM
         aux(k, 20) = Tg500_DM
         aux(k, 21) = n_e500
         aux(k, 22) = Ke500
         aux(k, 23) = Pe500
         aux(k, 24) = r200_DM
         aux(k, 25) = c200_DM
         aux(k, 26) = 0!Mg200_DM
         aux(k, 27) = MT200_DM
         aux(k, 28) = 0!fg200_DM
         aux(k, 29) = Tg200_DM
         aux(k, 30) = n_e200
         aux(k, 31) = Ke200
         aux(k, 32) = Pe200

         !i = 0
         !DO m = 33, 129, 8
         !   i = i + 1
         !   aux(k, m) = rgx(i)
         !   aux(k, m + 1) = M_DMx(i)
         !   aux(k, m + 2) = Mg_DMx(i)
         !   aux(k, m + 3) = fg_DMx(i)
         !   aux(k, m + 4) = Tgx(i)
         !   aux(k, m + 5) = n_ex(i)
         !   aux(k, m + 6) = Kex(i)
         !   aux(k, m + 7) = Pex(i)
         !END DO
         DEALLOCATE (rgx, Tgx, n_ex, Kex, Pex)
         DEALLOCATE (rhogasx, n_Hx, ne_nHx)
         DEALLOCATE (M_DMx, Mg_DMx, fg_DMx)

      END IF

      RETURN
   END SUBROUTINE MakeXrayGasDistributions

!================================================================================================

   FUNCTION DM_GNFWgasVol(radius, rs_DM, rp_GNFW, a_GNFW, b_GNFW, c_GNFW)

      IMPLICIT NONE

      REAL*8, PARAMETER   :: eps = 1.0d-4
      REAL*8               :: result
      REAL*8               :: radius, rs_DM, rp_GNFW, a_GNFW, b_GNFW, c_GNFW
      REAL*8               ::DM_GNFWgasVol
      CALL qtrap(DM_GNFWsphVolInt, 1.0d-4, radius, eps, result)
      DM_GNFWgasVol = result

   END FUNCTION DM_GNFWgasVol
!=========================================================================

   FUNCTION DM_GNFWsphVolInt(r)

      ! This is proportional to the density of the gas at radius r
      ! There is a constant coefficient not included here that makes
      ! it equal to pressure.
      !
      ! The function seems mostly used to calculate mass in a radius,
      ! by integrating over density
      !
      ! See Olamaie 2012, eq 6

      IMPLICIT NONE
      REAL*8               :: r
      REAL*8               ::DM_GNFWsphVolInt

      ! Undefined at zero?
      IF (r < r_min) THEN
         DM_GNFWsphVolInt = ((r_min*r_min*r_min)/((DLOG(1.0 + (r_min/rs_DM))) - (1.0/(1.0 + (rs_DM/r_min)))))* &
            ((r_min/rp_GNFW)**(-1.0*c_GNFW))* &
            ((1.0 + (r_min/rp_GNFW)**(a_GNFW))**(-1.0*(a_GNFW + b_GNFW - c_GNFW)/a_GNFW))* &
            ((b_GNFW*((r_min/rp_GNFW)**(a_GNFW))) + c_GNFW)

      ELSEIF (r > r_integration_max) THEN

         DM_GNFWsphVolInt = 0.d0

      ELSE

         DM_GNFWsphVolInt = ((r*r*r)/((DLOG(1.0 + (r/rs_DM))) - (1.0/(1.0 + (rs_DM/r)))))* &
            ((r/rp_GNFW)**(-1.0*c_GNFW))* &
            ((1.0 + (r/rp_GNFW)**(a_GNFW))**(-1.0*(a_GNFW + b_GNFW - c_GNFW)/a_GNFW))* &
            ((b_GNFW*((r/rp_GNFW)**(a_GNFW))) + c_GNFW)

      END IF

   END FUNCTION DM_GNFWsphVolInt
!================================================================================================
   FUNCTION GNFWmodel3D(r, rp_GNFW, a_GNFW, b_GNFW, c_GNFW, Pei)

      IMPLICIT NONE
      REAL*8, intent(in) :: r, rp_GNFW, a_GNFW, b_GNFW, c_GNFW, Pei
      REAL*8            :: GNFWmodel3D

      GNFWmodel3D = Pei/(((r/rp_GNFW)**(c_GNFW))* &
         ((1.d0 + ((r/rp_GNFW)**(a_GNFW)))** &
         ((b_GNFW - c_GNFW)/a_GNFW)))
   END FUNCTION GNFWmodel3D

!================================================================================================
   FUNCTION polyEstimateNumberDensity(radius)
      ! Estimate the number density of the gas by a method
      ! derived in Appendix C of Ghirardini2019

      implicit none
      real*8, intent(in) :: radius ! Mpc
      real*8 :: x, v, f0, result
      real*8 :: eps = 1.0d-4
      real*8 :: polyEstimateNumberDensity

      x = radius/r500_DM

      ! Invert our equation for temperature to find n_e at r500
      ! n0 = T0 ** (-1/Gamma0)
      ! Then f0 = f(x=1) = n0^Gamma0
      f0 = 1./T0_poly

      v = x**Gamma0*GammaR/(1.+Gamma0)

      CALL qtrap(polyhvIntegrand, 1.d0, x, eps, result)

      ! Technically it should be (f0*v0 + result)/v but v(x=1)=1
      ! While the paper has mu as the denominator of eq. C.12,
      ! running through the math shows it should be v instead.
      polyEstimateNumberDensity = ((f0 + result)/v)**(1/Gamma0)

   END FUNCTION polyEstimateNumberDensity
!================================================================================================
   FUNCTION polyhvIntegrand(zz)
      ! Integrand from eq C.12 of Ghirardini2019

      implicit none
      real*8 :: zz
      real*8 :: u, hx, v
      real*8 :: polyhvIntegrand

      v = zz**Gamma0*GammaR/(1 + Gamma0)
      u = (log(1.+c500_dm*zz) - c500_dm*zz/(1.+c500_dm*zz))/(log(1 + zz) - zz/(1 + zz))
      hx = -2*u*Gamma0/(T0_poly*zz**2*zz**GammaR*(Gamma0 + 1))

      polyhvIntegrand = v*hx

   END FUNCTION

!================================================================================================
   FUNCTION polyTemperature(radius, ne_r)
      ! Derived from eq. 6 of Ghirardini2019, neglecting the intrinsic scatter

      implicit none

      real*8, intent(in) :: radius, ne_r
      real*8 :: polyTemperature

      polyTemperature = exp(log(T0_poly) + Gamma0*log(ne_r*1.d-6) + GammaR*log(radius/r500_DM))*Tg500_DM

   END FUNCTION

!================================================================================================
   FUNCTION calcDMmass(rs_DM, rhos_DM, ri_DM)

      IMPLICIT NONE

      REAL*8              :: rs_DM, rhos_DM, ri_DM
      REAL*8              :: calcDMmass

      calcDMmass = 4.0*Pi*rhos_DM*rs_DM*rs_DM*rs_DM* &
         (DLOG(1.0 + (ri_DM/rs_DM)) - (1.0/(1.0 + (rs_DM/ri_DM))))

   END FUNCTION calcDMmass

!================================================================================================
   FUNCTION calcneDM(r, rs_DM, rp_GNFW, a_GNFW, b_GNFW, c_GNFW)

      IMPLICIT NONE

      REAL*8              :: r, rs_DM, rp_GNFW, a_GNFW, b_GNFW, c_GNFW
      REAL*8              ::calcneDM

      calcneDM = ((r)/((DLOG(1.0 + (r/rs_DM))) - (1.0/(1.0 + (rs_DM/r)))))* &
         ((r/rp_GNFW)**(-1.0*c_GNFW))* &
         ((1.0 + (r/rp_GNFW)**(a_GNFW))**(-1.0*(a_GNFW + b_GNFW - c_GNFW)/a_GNFW))* &
         ((b_GNFW*((r/rp_GNFW)**(a_GNFW))) + c_GNFW)

   END FUNCTION calcneDM

!================================================================
   FUNCTION EinastoDM_GNFWgasvol(radius, r_2_DM, alpha_Einasto, rp_GNFW, a_GNFW, &
      b_GNFW, c_GNFW)

      IMPLICIT NONE

      REAL*8, PARAMETER   :: eps = 1.0d-4
      REAL*8               :: result
      REAL*8               :: radius, r_2_DM, alpha_Einasto, rp_GNFW, a_GNFW, b_GNFW, c_GNFW
      REAL*8               :: EinastoDM_GNFWgasvol

      CALL qtrap(EinastoDM_GNFWsphVolInt, 1.0d-4, radius, eps, result)

      EinastoDM_GNFWgasVol = result

   END FUNCTION EinastoDM_GNFWgasvol
!=========================================================================

   FUNCTION EinastoDM_GNFWsphVolInt(r)

      IMPLICIT NONE
      REAL*8               :: r
      REAL*8               ::EinastoDM_GNFWsphVolInt

      IF (r < r_min) THEN
         EinastoDM_GNFWsphVolInt = ((r_min*r_min*r_min)* &
            ((b_GNFW*((r_min/rp_GNFW)**(a_GNFW))) + c_GNFW)* &
            ((r_min/rp_GNFW)**(-1.0*c_GNFW))* &
            (1.0d0 + (r_min/rp_GNFW)**(a_GNFW))** &
            (-1.0*(a_GNFW + b_GNFW - c_GNFW)/a_GNFW))/ &
            INCOG(Gamma_coeff1, ((2.d0/alpha_Einasto)*((r_min/r_2_DM)**alpha_Einasto)))
      ELSEIF (r > r_integration_max) THEN
         EinastoDM_GNFWsphVolInt = 0.d0
      ELSE
         EinastoDM_GNFWsphVolInt = ((r*r*r)* &
            ((b_GNFW*((r/rp_GNFW)**(a_GNFW))) + c_GNFW)* &
            ((r/rp_GNFW)**(-1.0*c_GNFW))* &
            (1.0d0 + (r/rp_GNFW)**(a_GNFW))** &
            (-1.0*(a_GNFW + b_GNFW - c_GNFW)/a_GNFW))/ &
            INCOG(Gamma_coeff1, ((2.d0/alpha_Einasto)*((r/r_2_DM)**alpha_Einasto)))
      END IF

   END FUNCTION EinastoDM_GNFWsphVolInt

!======================================================

   SUBROUTINE Xray_flux_coeff(Rhogas0, KB_Tx, n_e0, n_H0, ne_nH0, xrayE1, xrayE2, xrayNbin, xrayDeltaE, xrayFluxCoeff)

      IMPLICIT NONE
      INTEGER     :: i, k, flag, j, m, xrayNbin

      REAL*8       ::n_H0, Rhogas0, KB_Tx, xrayDeltaE
      REAL*8       ::Rhogas0_SI, el_den, n_e0
      REAL*8       ::ne_nH0, sum_xenuctot
      INTEGER, PARAMETER          :: NOEL = 15
      INTEGER, PARAMETER          :: NL_MAX = 5500
      REAL*8       ::xrayE1(1:xrayNbin), xrayE2(1:xrayNbin), xrayFluxCoeff(1:xrayNbin)

      REAL*8 xe, xzin, elx, flx
      DIMENSION xe(NOEL), xzin(NUMION), elx(NL_MAX), flx(NL_MAX)

      ! Initialise variables
      xrayde = 0.0
      DO i = 1, xrayNbin
         xrayFluxCoeff(i) = 0.
      END DO

      ! Convert gas density to SI units
      Rhogas0_SI = Rhogas0*m_sun/(Mpc2m*Mpc2m*Mpc2m)      !this is central gas density in kgm^-3

      ! Part of Eq 21 from the Bayes-X paper?
      sum_xenuctot = 0.
      DO i = 1, NOEL
         xe(i) = 10.**(abund(i) - 12.)
         sum_xenuctot = sum_xenuctot + xe(i)*nuc_tot(i)
      END DO

      n_H0 = Rhogas0_SI/(m_p_kg*sum_xenuctot)           !this is central hydrogen density in m^-3

      CALL MEKAL1(KB_Tx, xe, xzin, ne_nH0)

      el_den = ne_nH0*n_H0

      n_H0 = n_H0*1.d-6         ! central Hydrogen density in cm^-3
      n_e0 = el_den*1.d-6       ! central electron density in cm^-3
      el_den = n_e0              ! central electron density in cm^-3

      xfluxsec3 = (n_H0*n_H0)*(ne_nH0)/(SQRT(KB_Tx))

      CALL MEKAL4(xrayE1, xrayE2, xrayFluxcoeff, xrayNbin, KB_Tx, xzin)

      CALL MEKAL5(xrayE1(1), xrayE2(xrayNbin), KB_Tx, xzin, el_den, elx, flx, nlx)

      nl = nlx

      DO i = 1, xrayNbin
         xrayde = xrayE2(i) - xrayE1(i)
         DO WHILE (nl .GT. 0)
            IF (elx(nl) .GT. xrayE2(i)) GOTO 500
            xrayFluxCoeff(i) = xrayFluxCoeff(i) + flx(nl)/xrayde
            nl = nl - 1
         END DO
500   END DO

      DO j = 1, xrayNbin
         xrayFluxCoeff(j) = xfluxsec3*xrayFluxCoeff(j)*xrayDeltaE
      END DO

      RETURN
   END SUBROUTINE Xray_flux_coeff
!================================================================
   SUBROUTINE MEKAL1(T, Xe, Xzin, Ed)
      IMPLICIT NONE
      REAL*8 T, Ed
      REAL*8 alfa, fn, sion, so, Xe, Xzin
      INTEGER i, iel, j, k, nion, NIONMAX, NOEL

      PARAMETER(NOEL=15, NIONMAX=29)
      DIMENSION Xe(NOEL), Xzin(NUMION), alfa(NUMION), sion(NUMION), &
         fn(NIONMAX), nion(NOEL)
      DATA nion/2, 3, 7, 8, 9, 11, 12, 13, 14, 15, 17, 19, &
         21, 27, 29/

      CALL MEKAL2(sion, T)

      CALL MEKAL3(alfa, T)

      k = 2
      Xzin(1) = 0.0
      Xzin(2) = 1.0
      DO iel = 2, NOEL
         fn(1) = 1.0
         so = 1.0
         DO i = 2, nion(iel)
            fn(i) = fn(i - 1)*sion(k + i - 1)/alfa(k + i)
            IF (fn(i) .GT. 1.E25) THEN
               DO j = 1, i
                  fn(j) = fn(j)/1.E25
               END DO
               so = so/1.E25
            END IF
            so = so + fn(i)
         END DO
         DO i = 1, nion(iel)
            k = k + 1
            Xzin(k) = fn(i)/so
         END DO
      END DO
      Ed = 0.
      k = 0
      DO iel = 1, NOEL
         DO i = 1, nion(iel)
            k = k + 1
            IF (Xzin(k) .LT. 1E-5) THEN
               Xzin(k) = 0.
            ELSE
               Xzin(k) = Xzin(k)*Xe(iel)
               Ed = Ed + Xzin(k)*FLOAT(i - 1)
            END IF
         END DO
      END DO
      RETURN
   END SUBROUTINE MEKAL1
!================================================================

   SUBROUTINE MEKAL2(Sion, T)
      IMPLICIT NONE
      REAL*8 T, y, yay
      REAL*8 auto, f1, f2, f4, fb, fuy, g, &
         sa, sd, Sion, sk, t15, tsq, x, yfuy
      REAL*8 YLIM
      INTEGER ia, ind, j, k, n, nj
      INTEGER OKOK

      PARAMETER(YLIM=69.0)
      DIMENSION Sion(NUMION)

      t15 = T**1.5
      tsq = SQRT(T)
      ia = 0
      DO k = 1, NUMION
         Sion(k) = 0.
         ia = ia + 1
         nj = NINT(AIOn(ia))
         IF (nj .NE. 0) THEN
            sa = 0.
            DO j = 1, nj
               y = AIOn(ia + 1)/T
               x = 1./y
               IF (x .GT. 0.055) THEN
                  f4 = (-5.725E-4 + x*(0.01345 + x*(0.8691 + 0.03404*x))) &
                     /(1.+x*(2.197 + x*(0.2454 + 2.053E-3*x)))
               ELSE
                  f4 = .7699*x**1.9496
               END IF
               fuy = FMEKAL13(y)
               IF (y .LT. 10.0) THEN
                  fb = 1.+y - y*fuy*(2.+y)
               ELSE
                  CALL MEKAL12(y, fb)
               END IF
               IF (y .LT. YLIM) sa = sa + EXP(-y) &
                  *((AIOn(ia + 2)*(1.-y*fuy) + AIOn(ia + 3) &
                  *fb + AIOn(ia + 4)*fuy)/y + AIOn(ia + 5) &
                  *f4)
               ia = ia + 5
            END DO
            Sion(k) = sa/t15
         END IF
         ia = ia + 1
         ind = NINT(AIOn(ia))
         IF (ind .EQ. 1) THEN
         ELSEIF (ind .EQ. 3) THEN
            y = AIOn(ia + 1)/T
            IF (y .LT. YLIM) THEN
               fuy = FMEKAL13(y)
               yfuy = y*fuy
               sa = EXP(-y) &
                  *(AIOn(ia + 2)*fuy + AIOn(ia + 3)*(1.0 - yfuy) + AIOn(ia + 4) &
                  *yfuy + AIOn(ia + 5)*y*(1.0 - yfuy))
               Sion(k) = Sion(k) + sa/tsq
            END IF
            ia = ia + 5
         ELSEIF (ind .EQ. 4) THEN
            y = AIOn(ia + 1)/T
            IF (y .LT. YLIM) Sion(k) = Sion(k) + AIOn(ia + 2)*EXP(-y) &
               *(1.0 - y*FMEKAL13(y))/tsq
            ia = ia + 2
         ELSEIF (ind .EQ. 5) THEN
            sd = 0.
            DO n = 1, 18
               y = AIOn(ia + 1)/T
               IF (y .LT. YLIM) THEN
                  yay = y*AIOn(ia + 2)
                  fuy = FMEKAL13(yay)
                  sd = sd + EXP(-y) &
                     *(FMEKAL13(y)*AIOn(ia + 3) + AIOn(ia + 4) + y*(AIOn(ia + 5) &
                     *fuy + AIOn(ia + 6)*(1.-yay*fuy)))
               END IF
               ia = ia + 6
            END DO
            Sion(k) = Sion(k) + sd/tsq
         ELSEIF (ind .EQ. 6) THEN
            sk = 0.
            DO j = 1, 12
               y = AIOn(ia + 1)/T
               IF (y .LT. YLIM) THEN
                  f1 = EXP(-y)*(y + 1.)
                  y = AIOn(ia + 2)/T
                  IF (y .LT. YLIM) THEN
                     f2 = EXP(-y)*(y + 1.)
                     sk = sk + (f1 - f2)*AIOn(ia + 3)
                  END IF
               END IF
               ia = ia + 3
            END DO
            Sion(k) = Sion(k) + sk*tsq
            GOTO 100
         ELSEIF (ind .EQ. 7) THEN
            GOTO 100
         ELSE
            y = AIOn(ia + 1)/T
            IF (y .LT. YLIM) THEN
               auto = AIOn(ia + 2)*EXP(-y)*(1.+AIOn(ia + 3)*FMEKAL13(y))/tsq
               Sion(k) = Sion(k) + auto
            END IF
            ia = ia + 3
         END IF
         GOTO 200
100      y = AIOn(ia + 1)/T
         g = FMEKAL10(y, AIOn(ia + 2), AIOn(ia + 3), AIOn(ia + 4), AIOn(ia + 5), &
            AIOn(ia + 6))
         IF (y .LT. YLIM) Sion(k) = Sion(k) + g*EXP(-y)/tsq
         ia = ia + 6
200      OKOK = 0
      END DO
      IF (ia .GT. NAMAX .and. myID == 0) WRITE (*, *) 'INCOMPATIBLE DATASET'
      RETURN
   END SUBROUTINE MEKAL2
!===================================================================================
   SUBROUTINE MEKAL3(Alfa, T)
      IMPLICIT NONE
      REAL*8 T
      REAL*8 Alfa, diel, t13, t15, t25, tlo, tsq, y, YLIM
      INTEGER ind, ip, k, l, m

      PARAMETER(YLIM=69.0)
      DIMENSION Alfa(NUMION)
      tlo = DLOG(T)
      tsq = SQRT(T)
      t13 = T**0.33333333
      t15 = T**1.5
      t25 = T**2.5
      DO k = 1, NUMION
         ind = NINT(AREc(1, k))
         IF (ind .GT. 0) THEN
            IF (ind .LE. 8) THEN
               Alfa(k) = AREc(2, k)/T**AREc(3, k)
            ELSE
               Alfa(k) = AREc(2, k)/T**(AREc(3, k) + AREc(4, k)*T)
            END IF
            IF (ind .EQ. 2) THEN
               IF (T .GT. AREc(4, k)) THEN
                  Alfa(k) = Alfa(k) + AREc(5, k)/T**AREc(6, k)
               ELSE
                  y = AREc(7, k)/T
                  IF (y .LT. YLIM) THEN
                     diel = AREc(8, k)*EXP(-y)/t15
                     y = AREc(9, k)/T
                     IF (y .LT. YLIM) &
                        diel = diel*(1.+AREc(10, k)*EXP(-y))
                     Alfa(k) = Alfa(k) + diel
                  END IF
               END IF
            ELSEIF (ind .EQ. 3) THEN
               y = AREc(4, k)/T
               IF (y .LT. YLIM) THEN
                  diel = AREc(5, k)*EXP(-y)/t15
                  y = AREc(6, k)/T
                  IF (y .LT. YLIM) diel = diel*(1.+AREc(7, k)*EXP(-y)) &
                     **AREc(8, k)
                  y = AREc(9, k)/T
                  IF (y .LT. YLIM) diel = diel*(1.+AREc(10, k)*EXP(-y)) &
                     **AREc(11, k)
                  Alfa(k) = Alfa(k) + diel
               END IF
            ELSEIF (ind .EQ. 4) THEN
               IF (T .GT. AREc(4, k)) THEN
                  y = AREc(5, k)/T
                  IF (y .LT. YLIM) THEN
                     diel = AREc(6, k)*EXP(-y)/t15
                     y = AREc(7, k)/T
                     IF (y .LT. YLIM) diel = diel*(1.+AREc(8, k)*EXP(-y))
                     Alfa(k) = Alfa(k) + diel
                  END IF
               ELSE
                  Alfa(k) = Alfa(k) + AREc(9, k)/T**AREc(10, k)
               END IF
            ELSEIF (ind .EQ. 5) THEN
               IF (T .GT. AREc(4, k)) THEN
                  y = AREc(5, k)/T
                  IF (y .LT. YLIM) THEN
                     diel = AREc(6, k)*EXP(-y)/t15
                     y = AREc(7, k)/T
                     IF (y .LT. YLIM) diel = diel*(1.+AREc(8, k)*EXP(-y))
                     Alfa(k) = Alfa(k) + diel
                  END IF
               ELSE
                  y = AREc(9, k)/T
                  IF (y .LT. YLIM) Alfa(k) = Alfa(k) + AREc(10, k) &
                     /t15*EXP(-y)
               END IF
            ELSEIF (ind .EQ. 6) THEN
               y = AREc(4, k)/T
               IF (y .LT. YLIM) THEN
                  diel = AREc(5, k)*EXP(-y)/t15
                  y = AREc(6, k)/T
                  IF (y .LT. YLIM) diel = diel*(1.+AREc(7, k)*EXP(-y))
                  Alfa(k) = Alfa(k) + diel
               END IF
               IF (T .LT. AREc(8, k)) THEN
                  y = AREc(9, k)/T
                  IF (y .LT. YLIM) Alfa(k) = Alfa(k) + AREc(10, k) &
                     *EXP(-y) &
                     /t25*(((T + AREc(11, k))*T + AREc(12, k))*T + AREc(13, k))
               END IF
            ELSEIF (ind .EQ. 7) THEN
               y = AREc(4, k)/T
               IF (y .LT. YLIM) THEN
                  diel = AREc(5, k)*EXP(-y)/t15
                  y = AREc(6, k)/T
                  IF (y .LT. YLIM) diel = diel*(1.+AREc(7, k)*EXP(-y))
                  Alfa(k) = Alfa(k) + diel
               END IF
            ELSEIF (ind .EQ. 8 .OR. ind .EQ. 9) THEN
               ip = 4
               IF (ind .EQ. 9) ip = ip + 1
               DO l = 1, 4
                  m = ip + 2*l - 2
                  IF (AREc(m + 1, k) .NE. 0.) THEN
                     y = AREc(m, k)/T
                     IF (y .LT. YLIM) Alfa(k) = Alfa(k) + AREc(m + 1, k) &
                        /t15*EXP(-y)
                  END IF
               END DO
            ELSE
               y = AREc(4, k)/T
               IF (y .LT. YLIM) Alfa(k) = Alfa(k) + AREc(5, k) &
                  /t15*EXP(-y)
            END IF
         ELSEIF (ind .EQ. 0) THEN
            Alfa(k) = (AREc(2, k) - AREc(3, k)*tlo + AREc(4, k)*t13)/tsq
         ELSE
            Alfa(k) = 0.
         END IF
      END DO
   END SUBROUTINE MEKAL3
!=======================================================
   SUBROUTINE MEKAL4(E1, E2, Flx, Nemx, T, Xzin)
      IMPLICIT NONE

      REAL*8           ::T
      REAL*8 a, a31, a32, bf, CMAX, dgfbk, dx, e2p, ej, &
         f2p, g, g2a, g2p, gfb, gfbk, gff
      REAL*8 gl, phi, psi, twz, u, ua, ul, UMAX, w, W3PI, &
         x, x50, xez, xzg, Xzin, y
      INTEGER i, i1, i2, iel, ifo, in, in6, inmax, ip, iu, j, &
         j1, j2, j3, ju, k, k1, n, n0, n2
      INTEGER NA, NBF, nel, Nemx, nk, NO2P, NOEL, NP
      REAL*8           ::E1(Nemx), E2(Nemx), Flx(Nemx)

      PARAMETER(NOEL=15, NO2P=2)
      PARAMETER(NP=51, NBF=15, NA=NUU + 5 + 6*NBF)
      PARAMETER(UMAX=69.0, CMAX=1.E-3, W3PI=0.551328895)
      DIMENSION Xzin(NUMION), a(NA, NUMION)
      DIMENSION g2a(NG), ua(NUU), bf(4, NBF), phi(NP, NO2P)
      DIMENSION nel(NOEL), e2p(NOEL, NO2P), f2p(NOEL, NO2P)
      DATA bf/6.9248, -2.1368, 1.4878, -0.5867, 7.0540, -2.3572, &
         2.4469, -1.1122, 7.0799, -2.3877, 2.5579, -1.1304, &
         7.0873, -2.3932, 2.5715, -1.1090, 7.0865, -2.3896, &
         2.5563, -1.0801, 7.0837, -2.3845, 2.5371, -1.0550, &
         7.0769, -2.3769, 2.5122, -1.0299, 7.0728, -2.3724, &
         2.4966, -1.0133, 7.0699, -2.3692, 2.4848, -1.0007, &
         7.0658, -2.3638, 2.4690, -0.9871, 7.0572, -2.3588, &
         2.4556, -0.9764, 7.0530, -2.3551, 2.4443, -0.9667, &
         7.0557, -2.3555, 2.4422, -0.9626, 7.0487, -2.3492, &
         2.4259, -0.9510, 7.0489, -2.3481, 2.4206, -0.9455/

      DATA nel/2, 5, 12, 20, 29, 40, 52, 65, 79, 94, 111, &
         130, 151, 178, 207/
      DATA e2p/0.0102, 0.0408, 0.368, 0.500, 0.654, 1.022, 1.236, &
         1.473, 1.728, 2.006, 2.621, 3.324, 4.105, 6.965, &
         8.051, 0.0000, 0.0200, 0.304, 0.426, 0.569, 0.915, &
         1.119, 1.343, 1.575, 1.853, 2.450, 3.124, 3.887, &
         6.666, 7.798/

      DATA f2p/15*0.415, 0.000, 0.274, 0.645, 0.671, 0.691, &
         0.719, 0.729, 0.737, 0.745, 0.751, 0.761, 0.768, &
         0.775, 0.787, 0.790/

      DATA phi/0.000, 0.426, 0.768, 1.049, 1.281, 1.476, 1.642, &
         1.784, 1.905, 2.010, 2.101, 2.181, 2.252, 2.314, &
         2.367, 2.412, 2.451, 2.484, 2.513, 2.538, 2.559, &
         2.576, 2.589, 2.597, 2.602, 2.603, 2.602, 2.597, &
         2.589, 2.576, 2.559, 2.538, 2.513, 2.484, 2.451, &
         2.412, 2.367, 2.314, 2.252, 2.181, 2.101, 2.010, &
         1.905, 1.784, 1.642, 1.476, 1.281, 1.049, 0.768, &
         0.426, 0.000, 0.000, 0.122, 0.323, 0.617, 0.920, &
         1.176, 1.417, 1.600, 1.772, 1.923, 2.065, 2.177, &
         2.282, 2.378, 2.462, 2.541, 2.604, 2.661, 2.708, &
         2.742, 2.777, 2.803, 2.820, 2.832, 2.850, 2.859, &
         2.850, 2.832, 2.820, 2.803, 2.777, 2.742, 2.708, &
         2.661, 2.604, 2.541, 2.462, 2.378, 2.282, 2.177, &
         2.065, 1.923, 1.772, 1.600, 1.417, 1.176, 0.920, &
         0.617, 0.323, 0.122, 0.000/

      DO i = 1, NG
         g2a(i) = -4.25 + FLOAT(i)*0.25
      END DO
      DO i = 1, NUU
         ua(i) = -4.25 + FLOAT(i)*0.25
      END DO
      iel = 1
      k = 0
      DO i = 1, NUMION
         k1 = k + 1
         IF (k1 .LT. NUMION) THEN
            DO j = 1, NA
               a(j, k1) = 0.
            END DO
         END IF
         ifo = 0
         IF (Xzin(i) .GT. 0.) THEN
            xez = Xzin(i)*P(1, i)
            IF (xez .GT. CMAX) THEN
               ifo = 1
               k = k + 1
               a(1, k) = xez
               gl = DLOG10(.0136*P(1, i)/T)
               i1 = INT(4.*(gl + 4.)) + 1
               IF (i1 .LT. 1) i1 = 1
               IF (i1 .GE. NG) i1 = NG - 1
               i2 = i1 + 1
               w = (g2a(i2) - gl)/(g2a(i1) - g2a(i2))
               DO ju = 1, NUU
                  a(ju + 1, k) = GA(i2, ju) + w*(GA(i2, ju) - GA(i1, ju))
               END DO
            END IF
            DO i2 = 1, NO2P
               IF (i .EQ. nel(iel) - i2) THEN
                  IF (i2 .NE. 1) THEN
                     g = 0.05
                  ELSEIF (iel .EQ. 1) THEN
                     g = FMEKAL10((e2p(iel, i2)/T) &
                        , 0.08d0, -0.16d0, 0.11d0, 0.0d0, 0.0d0)
                  ELSE
                     g = 0.055
                  END IF
                  twz = 164995.*f2p(iel, i2)*g*Xzin(i)/e2p(iel, i2)
                  IF (twz .GT. CMAX) THEN
                     u = e2p(iel, i2)/T
                     IF (u .LT. UMAX) THEN
                        IF (ifo .EQ. 0) THEN
                           ifo = 1
                           k = k + 1
                        END IF
                        a(27, k) = i2
                        a(28, k) = e2p(iel, i2)
                        a(29, k) = twz*EXP(-u)
                     END IF
                  END IF
               END IF
            END DO
            IF (i .EQ. nel(iel)) iel = iel + 1
            n0 = NINT(P(2, i))
            IF (n0 .GT. 0) THEN
               in = 0
               DO n = n0, NBF
                  n2 = n*n
                  IF (n .EQ. n0) THEN
                     a31 = bf(1, n)*P(3, i)*Xzin(i)/T
                     a32 = P(4, i)/T
                  ELSE
                     a31 = bf(1, n)*P(5, i)*Xzin(i)/T/FLOAT(n2*n)
                     a32 = P(6, i)/T/FLOAT(n2)
                  END IF
                  IF (a32 .LT. UMAX) THEN
                     IF (a31*EXP(a32) .GT. CMAX) THEN
                        IF (ifo .EQ. 0) THEN
                           ifo = 1
                           k = k + 1
                        END IF
                        in6 = 6*in
                        a(31 + in6, k) = a31
                        a(32 + in6, k) = a32
                        a(33 + in6, k) = a32*FLOAT(n2)
                        a(34 + in6, k) = bf(2, n)
                        a(35 + in6, k) = bf(3, n)
                        a(36 + in6, k) = bf(4, n)
                        in = in + 1
                     END IF
                  END IF
               END DO
            END IF
            IF (ifo .EQ. 1) a(30, k) = in
         END IF
      END DO
      nk = k
      DO j = 1, Nemx
         gff = 0.
         gfb = 0.
         g2p = 0.
         ej = 0.5*(E1(j) + E2(j))
         u = ej/T
         IF (u .LT. UMAX) THEN
            ul = DLOG10(u)
            iu = INT(4.*(ul + 4.)) + 1
            IF (iu .LT. 1) iu = 1
            IF (iu .GE. NUU) iu = NUU - 1
            j1 = iu
            j2 = j1 + 1
            j3 = j2 + 1
            w = (ua(j2) - ul)/(ua(j1) - ua(j2))
            DO k = 1, nk
               IF (a(1, k) .GT. 0.) THEN
                  IF (u .LT. 1.E-4) THEN
                     xzg = a(1, k)*(EXP(a(2, k)) - W3PI*DLOG(u/1.E-4))
                  ELSE
                     xzg = a(1, k)*EXP(a(j3, k) + w*(a(j3, k) - a(j2, k)))
                  END IF
                  gff = gff + xzg
               END IF
               gfbk = 0.
               in = 0
               inmax = NINT(a(30, k))
               dgfbk = 1.
               DO WHILE (in .LT. inmax .AND. dgfbk .GT. CMAX*gfbk)
                  in6 = 6*in
                  IF (a(31 + in6, k) .GT. 0) THEN
                     y = u - a(32 + in6, k)
                     IF (y .GE. 0.0 .AND. y .LT. UMAX) THEN
                        x = 1./(1.+SQRT(u/a(33 + in6, k)))
                        dgfbk = (((a(36 + in6, k)*x + a(35 + in6, k))*x + a(34 + in6 &
                           , k))*x + 1.)*x*a(31 + in6, k)*EXP(-y)
                        gfbk = gfbk + dgfbk
                     END IF
                  END IF
                  in = in + 1
               END DO
               gfb = gfb + gfbk
               IF (a(29, k) .GT. 0.) THEN
                  IF (ej .LT. a(28, k)) THEN
                     i2 = NINT(a(27, k))
                     x = ej/a(28, k)
                     x50 = 50.*x
                     ip = INT(x50)
                     dx = x50 - FLOAT(ip)
                     ip = ip + 1
                     psi = phi(ip, i2) + dx*(phi(ip + 1, i2) - phi(ip, i2))
                     g2p = g2p + a(29, k)*x*psi
                  END IF
               END IF
            END DO
            gff = gff*EXP(-u)
         END IF
         Flx(j) = (gff + gfb + g2p)/ej
      END DO
      RETURN
   END SUBROUTINE MEKAL4
!==============================================================
   SUBROUTINE MEKAL5(E1, E2, Xkt, Xzin, Elden, Elx, Flx, Nlx)

      IMPLICIT NONE

      REAL*8 a, dcor, el, Elx, eta, f, fg, fgc, &
         Flx, g, t, tau, xz, Xzin
      INTEGER i, lp3, Nlx, nzz

      REAL*8         ::E1, E2, xkt, Elden, y

      DIMENSION Xzin(*), Elx(*), Flx(*)
      INTEGER*4 iaux

      t = Xkt*1.16048E+7
      tau = Xkt*1.16048
      eta = Elden*1E-12
      Nlx = 0
      DO i = 1, NE
         el = RP(1, i)
         IF (el .LT. E1) RETURN
         IF (el .LT. E2) THEN
            nzz = LP(4, i)
            xz = Xzin(nzz)
            IF (xz .GT. 0.) THEN
               a = RP(3, i)
               y = el*a/Xkt
               IF (y .LT. 40.) THEN
                  g = FMEKAL10(y, RP(4, i), RP(5, i), RP(6, i), RP(7, i), RP(8, i))
                  f = RP(2, i)
                  iaux = NRL(i)
                  CALL MEKAL6(f, g, a, fg, el, y, t, Xzin, nzz, RP(1, i), LP(1, i), &
                     TRAns(i), Elden, iaux, IDNr, CDRcafe, NCAfe)
                  fgc = EXP(-y)*xz*1.646E+5
                  IF (eta .GT. 1E-10) THEN
                     IF (ABS(RP(11, i)) .GT. 1E-10 .OR. ABS(RP(14, i)) &
                        .GT. 1E-10) THEN
                        lp3 = LP(3, i)
                        CALL MEKAL9(dcor, eta, tau, RP(1, i), lp3)
                        fgc = fgc*dcor
                     END IF
                  END IF
                  Nlx = Nlx + 1
                  Elx(Nlx) = el
                  Flx(Nlx) = fg*fgc/el
               END IF
            END IF
         END IF
      END DO
      RETURN
   END SUBROUTINE MEKAL5
!=================================================================
   SUBROUTINE MEKAL6(F, G, A, Fg, El, Y, T, Xzin, Nzz, Rp, Lp, Trans, Elden, Lnum, &
      Idnr, Cdrcafe, Ncafe)
      IMPLICIT NONE

      REAL*8 A, b6, bdr, bii, brf, bri, c1, c1f, cdr, Cdrcafe
      REAL*8 cii, ciid, dr, drs, El, ex, ex1, ex2, F
      REAL*8 f5, f6, Fg, fn, fr, G
      REAL*8 rii, Rp, rr, T, t2, t5, x, xcdr, Xzin, xzn, xzn1
      REAL*8 yy, z, z12, z4
      INTEGER ica, idr, idrs, iel, iex, ife, ihel, ii, img
      INTEGER isat, itr, iz, izion, Lnum, Ncafe, Nzz
      REAL*8       ::Elden, Y
      INTEGER*2 Lp(*)
      INTEGER*4 Idnr(Ncafe)
      DIMENSION Xzin(*), Rp(*), Cdrcafe(Ncafe)
      DIMENSION izion(15), cdr(19), bdr(19), b6(19), ciid(19), bri(15)
      CHARACTER*1 aster, dum1str
      CHARACTER*8 Trans

      DATA img, ica, ife/8, 13, 14/
      DATA aster/'*'/
      DATA izion/1, 2, 6, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20, 26, 28/
      DATA cdr/21.84, 47.09, 87.63, 31.55, 76.39, 114.60, 21.60, &
         36.63, 55.61, 15.29, 8.74, 3.79, 0.284, 0.728, 21.84, &
         31.55, 76.39, 114.60, 47.09/
      DATA bdr/0.3, 0.14, 0.06, 0.3, 0.14, 0.06, 0.3, 0.14, &
         0.06, 0.27, 0.26, 0.24, 0.22, 0.21, 0.3, 0.3, 0.14, &
         0.06, 0.14/
      DATA b6/5E-6, 3E-5, 5E-5, 5E-6, 3E-5, 5E-5, 5E-6, 3E-5, &
         5E-5, 5E-6, 5E-6, 5E-6, 5E-6, 5E-6, 5E-6, 5E-6, &
         3E-5, 5E-5, 3E-5/
      DATA ciid/0.0, 0.0, 0.0, 0.106, 0.0, 0.0, 0.27, 0.0, 0.0, &
         0.28, 0.28, 0.46, 0.33, 0.3, 0.0, 0.106, 0.0, 0.0, &
         0.0/
      DATA bri/0.0, 0.0, 0.111, 0.225, 0.293, 0.338, 0.353, &
         0.370, 0.391, 0.425, 0.506, 0.595, 0.675, 0.821, 0.838/
      iel = Lp(1)
      iz = Lp(2) - 1
      z = FLOAT(izion(iel))
      z4 = z**4
      iex = Lp(5)
      idrs = Lp(6)
      idr = Lp(7)
      ii = Lp(8)
      itr = Lp(9)
      Fg = F*G
      t5 = T*1E-5
      dum1str = Trans(1:1)
      IF (iex .NE. 0 .AND. iex .NE. 2) THEN
         yy = (Y + 0.6)/(Y + 1.15)
         c1f = 1.+7.*yy*(1.065*(1 - bri(iel)) + 0.4*EXP(-0.21*Y))
         brf = FMEKAL7(z, t5, Elden, bri(iel))
         IF (iex .EQ. 3) THEN
            ex = (bri(iel) + c1f*(1.-brf)/yy/7.46)*1.065
         ELSE
            ex = c1f*brf
         END IF
      ELSE
         ex = 1.
      END IF
      ex = ex*Fg
      IF (idrs .NE. 0) THEN
         iz = iz + 1
         z12 = (iz + 1)**2
         x = (El*A)/(iz + 1)
         c1 = EXP(Y/(1.+z12/iz**3/.019))
         c1 = 8538.1*F*c1/((x*81.160 + 7.7235)*x + 1.)/T
         drs = c1*(A*El)**1.5*SQRT(iz/(iz*iz + 13.4))*z12
      ELSE
         drs = 0.
      END IF
      xzn1 = Xzin(Nzz + 1)
      xzn = Xzin(Nzz)
      dr = 0.
      rr = 0.
      IF (xzn1 .GT. 0) THEN
         IF (idr .NE. 0) THEN
            IF (idr .EQ. 1) THEN
               IF (iel .EQ. img .OR. iel .EQ. ica .OR. iel .EQ. ife) THEN
                  IF (dum1str .EQ. aster) THEN
                     DO isat = 1, Ncafe
                        IF (Idnr(isat) .EQ. Lnum) GOTO 4
                     END DO
4                    xcdr = Cdrcafe(isat)
                     GOTO 10
                  END IF
               END IF
               xcdr = cdr(itr)
10             dr = xcdr/z4*EXP(bdr(itr)*Y)/(b6(itr) + 1./(z - 1.)**4)
               rr = Rp(9)*T**Rp(10)
            ELSEIF (idr .EQ. 2) THEN
               dr = 11.*EXP(0.284*Y)/(1.+z4*6E-6) + 27.*EXP(0.1*Y)/(1.+z4*3E-5)
               rr = Rp(9)*T**Rp(10)
            ELSE
               ex1 = EXP(0.284*Y)
               ex2 = EXP(0.1*Y)
               f5 = 2.3*ex1 + 215.*ex2/(1.+z4*1.9E-4)
               f6 = 16.*ex1/(1.+z4*7E-5) + 90.*ex2/(1.+z4*8E-5)
               brf = FMEKAL7(z, t5, Elden, bri(iel))
               dr = FMEKAL8(bri(iel), brf, f5, f6, idr + 2)
               f5 = (iz + 1.)**2.62/T**0.31*18.72E-4
               f6 = (iz + 1.)**2.08/T**0.04*.575E-4
               rr = FMEKAL8(bri(iel), brf, f5, f6, idr + 2)
            END IF
            fr = xzn1/xzn*El/12.3985
            rr = fr*EXP(Y)*rr
            dr = fr*0.47*z4/T*dr
         ELSE
            rr = xzn1/xzn*EXP(Y)*El/12.3985*Rp(9)*T**Rp(10)
         END IF
      END IF
      rii = 0.
      IF (ii .NE. 0) THEN
         xzn1 = Xzin(Nzz - 1)
         xzn = Xzin(Nzz)
         IF (xzn1 .GT. 0) THEN
            IF (ii .EQ. 2) THEN
               IF (itr .NE. 0) THEN
                  cii = ciid(itr)
                  IF (itr .GE. 7 .AND. itr .LE. 14) THEN
                     bii = 1.242 + 0.009*(z - 26)
                  ELSE
                     bii = 1.296
                  END IF
               ELSE
                  bii = 1.164 + 0.008*(z - 26)
                  cii = 0.3
               END IF
               t2 = cii*FMEKAL11(bii*Y)*1.231/bii/(1.+2.E5*z**(-3.5))
            ELSE
               yy = (z - 2.)**4.3*6E12
               fn = yy + 1.33333333*Elden
               f5 = Elden/fn
               f6 = (yy + Elden/3.)/fn
               brf = FMEKAL7(z, t5, Elden, bri(iel))
               ihel = 5
               IF (ii .EQ. 1) ihel = 6
               t2 = 0.22*FMEKAL11(1.33*Y)*FMEKAL8(bri(iel), brf, f5, f6, ihel)
            END IF
            rii = xzn1/xzn*EXP(Y)*t2
         END IF
      END IF
      Fg = ex + drs + dr + rii + rr
      RETURN
   END SUBROUTINE MEKAL6
!====================================================================================

   FUNCTION FMEKAL7(Z, T5, Elden, Bri)
      IMPLICIT NONE

      REAL*8 betac, Bri, ex, FMEKAL7, T5, Z
      REAL*8       ::Elden
      ex = -0.73/Z**0.415
      betac = 330./Z**14.57*(T5/Z**2)**ex
      FMEKAL7 = 1./(1.+(betac*Elden)*Bri)
      RETURN
   END FUNCTION FMEKAL7

!================================================
   FUNCTION FMEKAL8(Bri, Brf, F5, F6, Itr)
      IMPLICIT NONE

      REAL*8 Brf, Bri, F5, F6, FMEKAL8
      INTEGER Itr

      IF (Itr .NE. 6) THEN
         FMEKAL8 = F5*Bri + (F6 + (1.-Bri)*F5)*(1.-Brf)
      ELSE
         FMEKAL8 = (F6 + (1.-Bri)*F5)*Brf
      END IF
      RETURN
   END FUNCTION FMEKAL8

!===========================================================
   SUBROUTINE MEKAL9(Dcor, Eta, Tau, Rp, Lp3)
      IMPLICIT NONE

      REAL*8 Dcor, Eta, ex, Rp, Tau
      INTEGER Lp3

      DIMENSION Rp(*)
      IF (Lp3 .EQ. 5) THEN
         Dcor = 1.+Rp(11)/(Eta**0.03 + 120./Tau**0.3/Eta)
      ELSEIF (ABS(Rp(14)) .GT. 1E-10) THEN
         Dcor = 1.+Rp(14)/(1.+3.5/Eta**0.75)
      ELSE
         ex = -Rp(12)*Eta**Rp(13)
         Dcor = 1.+Rp(11)*(1.-EXP(ex))
      END IF
      RETURN
   END SUBROUTINE MEKAL9


   FUNCTION FMEKAL10(Y, A, B, C, D, E)
      IMPLICIT NONE

      REAL*8 A, B, C, D, E, e1, f2, f3, f4, f5
      REAL*8 FMEKAL10
      REAL*8 Y
      IF (Y .LT. 1.0) THEN
         IF (Y .LT. 0.3) THEN
            e1 = -.57731566 + Y*(.99999193 + Y*(-.24991055 + Y*.05519968)) &
               - DLOG(Y)
            f5 = e1*DEXP(Y)
            f2 = Y*f5
         ELSE
            f2 = (.06225196 + Y*(1.646421 + Y*1.040425)) &
               /(.85539 + Y*(2.754082 + Y))
            f5 = f2/Y
         END IF
         f3 = Y*(1.-f2)
         f4 = Y*(1.-f3)
      ELSE
         IF (Y .LT. 10.) THEN
            f2 = (.250621 + Y*(2.334733 + Y))/(1.681534 + Y*(3.330657 + Y))
            f3 = Y*(1.-f2)
            f4 = Y*(1.-f3)
         ELSE
            f4 = (.37559 + Y*(-.6993774 + Y*2.))/(-2.520804 + Y*(2.617596 + Y))
            f3 = 1.-f4/Y
            f2 = 1.-f3/Y
         END IF
         f5 = f2/Y
      END IF
      FMEKAL10 = A + B*f2 + C*f3 + D*f4 + E*f5
      RETURN
   END FUNCTION FMEKAL10


   FUNCTION FMEKAL11(X)
      IMPLICIT NONE

      REAL*8 a, FMEKAL11
      REAL*8 X

      DIMENSION a(6)
      DATA a/-0.57721566, 0.99999193, -0.24991055, 0.05519968, -0.00976004, 0.00107857/
      IF (X .LE. 1.) THEN
         FMEKAL11 = -DLOG(X) + ((((a(6)*X + a(5))*X + a(4))*X + a(3))*X + a(2))*X + a(1)
      ELSE
         FMEKAL11 = (.250621/X + X + 2.334733)*EXP(-X)/(X*X + 3.330657*X + 1.681534)
      END IF
      RETURN
   END FUNCTION FMEKAL11


   SUBROUTINE MEKAL12(Y, Fb)
      IMPLICIT NONE

      REAL*8 c, Fb, ff, sb, t
      REAL*8       ::Y
      INTEGER i

      DIMENSION t(16), c(16)
      DATA t/.8764941047E-1, .4626963289, 1.141057774, 2.129283645, &
         3.437086633, 5.078018614, 7.070338535, 9.438314336, &
         12.21422336, 15.44152736, 19.18015685, 23.51590569, &
         28.57872974, 34.58339870, 41.94045264, 51.70116033/
      DATA c/.2061517149, .3310578549, .2657957776, .1362969342, &
         .4732892869E-1, .1129990008E-1, .1849070943E-2, &
         .2042719153E-3, .1484458687E-4, .6828319330E-6, &
         .1881024841E-7, .2862350242E-9, .2127079033E-11, &
         .6297967002E-14, .50504737E-17, .416146237E-21/
      sb = 0.
      DO i = 1, 16
         ff = t(i)/(t(i) + Y)
         sb = sb + ff**2*c(i)
      END DO
      Fb = sb
      RETURN
   END SUBROUTINE MEKAL12

   FUNCTION FMEKAL13(U)
      IMPLICIT NONE

      REAL*8 FMEKAL13
      REAL*8 U

      IF (U .LT. 1.0) THEN
         FMEKAL13 = DLOG(1.+1./U) - (0.36 + 0.03/SQRT(U + 0.01))/(1.+U)**2
      ELSE
         FMEKAL13 = DLOG(1.+1./U) - (0.36 + 0.03*SQRT(U + 0.01))/(1.+U)**2
      END IF
      RETURN
   END FUNCTION FMEKAL13

   SUBROUTINE xsphab(ear, ne, N_H_col, photar)

      INTEGER ne
      ! REAL ear(0:ne), param(1), photar(ne)
      REAL param(1), photar(ne)
      REAL N_H_col
      REAL*8 ear(0:ne)

      ! "Wisconsin" absorption using the cross-sections of Balucinska-Church
      ! and McCammon, 1992, ApJ 400, 599. The elemental abundances are those
      ! set by the abund command.

      ! Arguments :
      !     ear       r        i: the energy ranges on which to calculate the model
      !     ne        i        i: the number of energy ranges
      !     param     r        i: the H column density in cm^-2
      !     photar    r        r: fractional transmission

      INTEGER NPARM
      PARAMETER(NPARM=19)

      REAL vparam(NPARM)

      INTEGER i

      param(1) = N_H_col*1.e-22 ! 10^22 cm^-2
      vparam(1) = param(1)
      DO i = 2, NPARM - 1
         vparam(i) = 1. ! relative values?
      END DO
      vparam(NPARM) = 0.

      CALL xszvph(ear, ne, vparam, photar)

      RETURN
   END SUBROUTINE xsphab

   SUBROUTINE xszvph(ear, ne, param, photar)

      INTEGER ne
      ! REAL ear(0:ne), param(19), photar(ne), photer(ne)
      REAL param(19), photar(ne), photer(ne)
      REAL*8 ear(0:ne)

      ! calculates the photoelectric absorption using the cross-sections of
      ! Balucinska-Church and McCammon, 1992, ApJ 400, 599. The elemental
      ! abundances are specified relative to the ratios set using the abund
      ! command.
      ! Parameters :
      !       1     H column in 10^22 cm^-2
      !       2-18  Relative abundances of
      !               He,C,N,O,Ne,Na,Mg,Al,Si,S,Cl,Ar,Ca,Cr,Fe,Co,Ni
      !       19    Redshift

      ! Arguments :
      !      ear     r        i: energy ranges
      !      ne      i        i: number of energies
      !      param   r        i: model parameters
      !      ifl     i        i: file number
      !      photar  r        o: transmitted fraction

      INTEGER NPARM
      PARAMETER(NPARM=19)

      REAL pparam(NPARM)

      INTEGER i

      pparam(1) = param(1)
      DO i = 2, NPARM - 1
         pparam(i) = param(i)*param(1)
      END DO
      pparam(NPARM) = param(NPARM)

      CALL xszvab(ear, ne, pparam, photar)

      RETURN
   END SUBROUTINE xszvph

   SUBROUTINE xszvab(ear, ne, param, photar)

      INTEGER ne
      ! REAL ear(0:ne), param(19), photar(ne)
      REAL param(19), photar(ne)
      REAL*8 ear(0:ne)

      !  works out the redshifted transmission for material whose abundances
      !  are specified for the following 18 elements :
      !    1   hydrogen    (1)
      !    2   helium      (2)
      !    3   carbon      (6)
      !    4   nitrogen    (7)
      !    5   oxygen      (8)
      !    6   neon        (10)
      !    7   sodium      (11)
      !    8   magnesium   (12)
      !    9   aluminium   (13)
      !   10   silicon     (14)
      !   11   sulphur     (16)
      !   12   chlorine    (17)
      !   13   argon       (18)
      !   14   calcium     (20)
      !   15   chromium    (24)
      !   16   iron        (26)
      !   17   cobalt      (27)
      !   18   nickel      (28)
      !  The parameters are the column densities of the 18 elements in units of
      !  the column of each element in a solar abundance column of equivalent
      !  hydrogen column density of 1e22 /cm/cm. Parameter 19 is the redshift.

      ! Arguments :
      !      ear     r        i: energy ranges
      !      ne      i        i: number of energies
      !      param   r        i: model parameters
      !      ifl     i        i: file number
      !      photar  r        o: transmitted fraction

      INTEGER NELTS
      PARAMETER(NELTS=18)
      REAL LYLIMIT
      PARAMETER(LYLIMIT=0.0135984)

      INTEGER atomic(NELTS)
      INTEGER ie, status, i

      REAL column(NELTS)
      REAL zfac, elow, ehi, xsect, pret

      CHARACTER*2 celts(NELTS)

      ! External references to cross section function

      DATA atomic/1, 2, 6, 7, 8, 10, 11, 12, &
         13, 14, 16, 17, 18, 20, 24, 26, &
         27, 28/

      DATA celts/'H ', 'He', 'C ', 'N ', 'O ', 'Ne', 'Na', 'Mg', &
         'Al', 'Si', 'S ', 'Cl', 'Ar', 'Ca', 'Cr', 'Fe', &
         'Co', 'Ni'/

      status = 0

      ! Set the element columns

      DO i = 1, NELTS
         column(i) = param(i)*1.e22*fgabnd(celts(i))
      END DO

      zfac = 1.0 + param(19) ! TODO: This isn't factoring in redshift correctlyS I think as param(19) is always 0
      elow = ear(0)*zfac

      DO ie = 1, ne

         ehi = ear(ie)*zfac
         xsect = 0.

         ! Set transmission to unity above 100 keV and below the Lyman limit

         IF (elow .GT. 100. .OR. elow .LT. LYLIMIT) THEN

            photar(ie) = 1.0

            ! else calculate the cross-section

         ELSE

            DO i = 1, NELTS

               IF (fgxsct() .EQ. 'bcmc') THEN
                  pret = photo(elow, ehi, atomic(i), 3, status)
               ELSEIF (fgxsct() .EQ. 'obcm') THEN
                  pret = photo(elow, ehi, atomic(i), 2, status)
               END IF
               IF (column(i) .GT. 0.) THEN
                  IF (pret .LT. 1e30/column(i)) THEN
                     xsect = xsect + pret*column(i)
                  ELSE
                     xsect = 1e30
                  END IF
               END IF

            END DO

            photar(ie) = EXP(-xsect)

         END IF

         elow = ehi

      END DO

      RETURN
   END SUBROUTINE xszvab

   FUNCTION fgabnd(el)
      REAL fgabnd
      CHARACTER*2 el
      IF (el == 'H ') THEN
         fgabnd = 1.00e+00
      ELSEIF (el == 'He') THEN
         fgabnd = 9.77e-02
      ELSEIF (el == 'Li') THEN
         fgabnd = 1.45e-11
      ELSEIF (el == 'Be') THEN
         fgabnd = 1.41e-11
      ELSEIF (el == 'B ') THEN
         fgabnd = 3.98e-10
      ELSEIF (el == 'C ') THEN
         fgabnd = 3.63e-04
      ELSEIF (el == 'N ') THEN
         fgabnd = 1.12e-04
      ELSEIF (el == 'O ') THEN
         fgabnd = 8.51e-04
      ELSEIF (el == 'F ') THEN
         fgabnd = 3.63e-08
      ELSEIF (el == 'Ne') THEN
         fgabnd = 1.23e-04
      ELSEIF (el == 'Na') THEN
         fgabnd = 2.14e-06
      ELSEIF (el == 'Mg') THEN
         fgabnd = 3.80e-05
      ELSEIF (el == 'Al') THEN
         fgabnd = 2.95e-06
      ELSEIF (el == 'Si') THEN
         fgabnd = 3.55e-05
      ELSEIF (el == 'P ') THEN
         fgabnd = 2.82e-07
      ELSEIF (el == 'S ') THEN
         fgabnd = 1.62e-05
      ELSEIF (el == 'Cl') THEN
         fgabnd = 3.16e-07
      ELSEIF (el == 'Ar') THEN
         fgabnd = 3.63e-06
      ELSEIF (el == 'K ') THEN
         fgabnd = 1.32e-07
      ELSEIF (el == 'Ca') THEN
         fgabnd = 2.29e-06
      ELSEIF (el == 'Sc') THEN
         fgabnd = 1.26e-09
      ELSEIF (el == 'Ti') THEN
         fgabnd = 9.77e-08
      ELSEIF (el == 'V ') THEN
         fgabnd = 1.00e-08
      ELSEIF (el == 'Cr') THEN
         fgabnd = 4.68e-07
      ELSEIF (el == 'Mn') THEN
         fgabnd = 2.45e-07
      ELSEIF (el == 'Fe') THEN
         fgabnd = 4.68e-05
      ELSEIF (el == 'Co') THEN
         fgabnd = 8.32e-08
      ELSEIF (el == 'Ni') THEN
         fgabnd = 1.78e-06
      ELSEIF (el == 'Cu') THEN
         fgabnd = 1.62e-08
      ELSEIF (el == 'Zn') THEN
         fgabnd = 3.98e-08
      END IF
      RETURN
   END FUNCTION fgabnd

   FUNCTION fgxsct()
      CHARACTER fgxsct*4
      fgxsct = 'bcmc'
      RETURN
   END FUNCTION fgxsct

   FUNCTION photo(keV1, keV2, Z, versn, status)
      !  Routine to return the photoelectric absorption cross section in cm**2
      !  at energy keV for element Z.

      REAL photo, keV1, keV2
      INTEGER Z, versn, status

      !  Cross-section data from Henke etal, Atomic Data and Nuclear Data Tables
      !  vol 27, no 1, 1982. Fits mainly by Monika Balucinska-Church and Dan McCammon
      !  "Photoelectric Absorption Cross Sections with Variable Abunances"
      !  Ap.J. 400, 699 (1992)

      !  kaa   4/16/92
      !       11/26/99     added version number so that new and old versions of
      !                    the He cross-section can be accomodated.

      !  Arguments :
      ! keV1 r i: Lower energy of bin in keV.
      ! keV2 r i: Upper energy of bin in keV.
      !       Z       i i: Atomic number of element
      !       versn   i       i: 2 == old Marr & West He x-section
      !                          3 == new Yan et al. x-section
      ! status i r: 0 = OK
      !      1 = No data for this element
      ! photo r r: Cross-section in cm**2

      ! Tabulated fit coefficients - coeff(i,j,k) is the ith coefficient for the
      ! jth energy range for the kth element. i goes from 1 to 9 and j from 1 to
      ! 4.

      INTEGER MAXEDG, MAXCOF, MAXZ
      PARAMETER(MAXEDG=3, MAXCOF=9, MAXZ=79)

      DOUBLE PRECISION coeff(MAXCOF, MAXEDG + 1, MAXZ)
      DOUBLE PRECISION edge(MAXEDG, MAXZ)
      DOUBLE PRECISION Elog, E1, E2, X, Xt

      REAL f(MAXZ)

      INTEGER ncoefs(MAXEDG + 1, MAXZ)
      INTEGER nedges(MAXZ)
      INTEGER i, j, i1, i2, k

      ! ***** Coefficients for polynomial fits to cross-sections *****

      ! Hydrogen

      DATA((coeff(i, j, 1), i=1, 9), j=1, 4)/ &
         21.46941, 0.9398479, -0.1492932, 5.4634294d-3, &
         32*0./

      ! Helium

      DATA((coeff(i, j, 2), i=1, 9), j=1, 4)/ &
         14.61546, 4.682793, -0.7323856, 4.6526663d-2, -1.1172282d-3, &
         31*0./

      ! Lithium

      DATA((coeff(i, j, 3), i=1, 9), j=1, 4)/36*0./

      ! Beryllium

      DATA((coeff(i, j, 4), i=1, 9), j=1, 4)/36*0./

      ! Boron

      DATA((coeff(i, j, 5), i=1, 9), j=1, 4)/36*0./

      ! Carbon

      DATA((coeff(i, j, 6), i=1, 9), j=1, 4)/ &
         8.74161, 7.13348, -1.14604, 0.0677044, 5*0., &
         3.81334, 8.93626, -1.06905, 0.0422195, &
         23*0./

      ! Nitrogen

      DATA((coeff(i, j, 7), i=1, 9), j=1, 4)/ &
         9.24058, 7.02985, -1.08849, 0.0611007, 5*0., &
         -13.0353, 15.4851, -1.89502, 0.0769412, &
         23*0./

      ! Oxygen

      DATA((coeff(i, j, 8), i=1, 9), j=1, 4)/ &
         2.57264, 10.9321, -1.79383, 0.102619, 5*0., &
         16.53869, 3.6428144, -0.3177744, 7.9471897d-3, &
         23*0./

      ! Fluorine

      DATA((coeff(i, j, 9), i=1, 9), j=1, 4)/36*0./

      ! Neon

      DATA((coeff(i, j, 10), i=1, 9), j=1, 4)/ &
         -3.04041, 13.0071, -1.93205, 0.0977639, 5*0., &
         17.6007, 3.29278, -0.263065, 5.68290d-3, &
         23*0./

      ! Sodium

      DATA((coeff(i, j, 11), i=1, 9), j=1, 4)/ &
         -2737.598, 2801.704, -1009.892, 87.16455, 43.20644, &
         -15.27259, 2.180531, -0.1526546, 4.3137977d-3, &
         1.534019, 9.261744, -0.9914126, 3.5278253d-2, &
         23*0./

      ! Magnesium

      DATA((coeff(i, j, 12), i=1, 9), j=1, 4)/ &
         7.107172, 3.7359418, 7*0., &
         -81.32915, 65.2775, -15.00826, 1.558686, -6.1339621d-2, 4*0., &
         -9.161526, 13.07448, -1.435878, 5.2728362d-2, &
         14*0./

      ! Aluminium

      DATA((coeff(i, j, 13), i=1, 9), j=1, 4)/ &
         26.90487, -6.135221, 1.175546, 6*0., &
         -38.1232, 29.5161, -4.45416, 0.226204, 5*0., &
         14.6897, 4.22743, -0.344185, 8.18542d-3, &
         14*0./

      ! Silicon

      DATA((coeff(i, j, 14), i=1, 9), j=1, 4)/ &
         -3.066295, 10.006248, -0.9627411, 6*0., &
         -182.7217, 128.061, -29.47269, 3.03284, -0.1173096, 4*0., &
         -33.39074, 21.42992, -2.385117, 8.887583d-2, &
         14*0./

      ! Phosphorus

      DATA((coeff(i, j, 15), i=1, 9), j=1, 4)/36*0./

      ! Sulphur

      DATA((coeff(i, j, 16), i=1, 9), j=1, 4)/ &
         598.2911, -675.2265, 308.1133, -68.99324, 7.62458, &
         -0.3335031, 3*0, &
         3994.831, -3690.886, 1417.287, -287.9909, 32.70061, &
         -1.968987, 4.9149349d-2, 2*0., &
         -22.49628, 17.24599, -1.848444, 6.6506132d-2, &
         14*0./

      ! Chlorine

      DATA((coeff(i, j, 17), i=1, 9), j=1, 4)/ &
         6253.247, -8222.248, 4491.675, -1302.145, 211.4881, &
         -18.25547, 0.6545154, 2*0., &
         -233.0502, 146.9776, -31.12463, 2.938618, -0.104096, &
         4*0., &
         -23.74675, 17.50997, -1.857953, 6.6208832d-2, &
         14*0./

      ! Argon

      DATA((coeff(i, j, 18), i=1, 9), j=1, 4)/ &
         -330.3509, 270.7433, -78.90498, 10.35983, -0.5140201, &
         4*0., &
         -5.71870, 8.85812, -0.307357, 0.00169351, -0.0138134, &
         0.00120451, 3*0, &
         19.1905, 2.74276, -0.164603, 0.00165895, &
         14*0./

      ! Potassium

      DATA((coeff(i, j, 19), i=1, 9), j=1, 4)/36*0./

      ! Calcium

      DATA((coeff(i, j, 20), i=1, 9), j=1, 4)/ &
         -873.972, 868.5231, -339.678, 66.83369, -6.590398, &
         0.2601044, 3*0, &
         -3449.707, 2436.409, -682.0668, 95.3563, -6.655018, &
         0.1854492, 3*0, &
         18.89376, 2.709646, -0.1377201, &
         15*0./

      ! Scandium

      DATA((coeff(i, j, 21), i=1, 9), j=1, 4)/36*0./

      ! Titanium

      DATA((coeff(i, j, 22), i=1, 9), j=1, 4)/36*0./

      ! Vanadium

      DATA((coeff(i, j, 23), i=1, 9), j=1, 4)/36*0./

      ! Chromium

      DATA((coeff(i, j, 24), i=1, 9), j=1, 4)/ &
         -0.4919405, 15.66939, -5.199775, 1.086566, -0.1196001, &
         5.2152011d-3, 3*0, &
         27.29282, 0.2966640, 7*0., &
         -15.2525, 16.23729, -1.966778, 8.062207d-2, 5*0., &
         8.307041, 5.008987, -0.2580816, 6*0./

      ! Manganese

      DATA((coeff(i, j, 25), i=1, 9), j=1, 4)/36*0./

      ! Iron

      DATA((coeff(i, j, 26), i=1, 9), j=1, 4)/ &
         -15.07332, 21.94335, -4.862457, 0.5573765, -3.0065542d-2, &
         4.9834867d-4, 3*0, &
         -253.0979, 138.4238, -25.47119, 2.08867, -6.4264648d-2, &
         4*0., &
         -1.037655, 7.022304, -0.3638919, &
         15*0./

      ! Cobalt

      DATA((coeff(i, j, 27), i=1, 9), j=1, 4)/ &
         9.171919, 3.5721176, 7*0., &
         -6.910097, 13.58385, -1.873453, 9.1612935d-2, 5*0., &
         13.96877, 2.128918, 0.1149042, 4.9106661d-02, -1.4725224d-02, &
         8.3086651d-04, 3*0., &
         28.72910, 0.4456830, 7*0./

      ! Nickel

      DATA((coeff(i, j, 28), i=1, 9), j=1, 4)/ &
         -7.919931, 14.06475, -1.935318, 9.3929626d-2, 5*0., &
         3.71129, 8.45098, -0.896656, 0.0324889, 5*0., &
         28.4989, 0.485797, &
         16*0./

      ! Gold   !! Important : this is only valid above 220 eV !!

      DATA((coeff(i, j, 79), i=1, 9), j=1, 4)/ &
         -27.40668, 18.11780, -1.869548, 6.2878355D-02, 5*0., &
         26.50165, 0.6361220, 7*0., &
         -33.83069, 19.80218, -1.989242, 6.7341216D-02, &
         14*0./

      ! Number of coefficients (-1 => range does not exist)
      ! H - Ne

      DATA((ncoefs(i, j), i=1, 4), j=1, 10)/ &
         4, -1, -1, -1, &
         5, -1, -1, -1, &
         -1, -1, -1, -1, &
         -1, -1, -1, -1, &
         -1, -1, -1, -1, &
         4, 4, -1, -1, &
         4, 4, -1, -1, &
         4, 4, -1, -1, &
         -1, -1, -1, -1, &
         4, 4, -1, -1/

      ! Na-Ar

      DATA((ncoefs(i, j), i=1, 4), j=11, 18)/ &
         9, 4, -1, -1, &
         2, 5, 4, -1, &
         3, 4, 4, -1, &
         3, 5, 4, -1, &
         -1, -1, -1, -1, &
         6, 7, 4, -1, &
         7, 5, 5, -1, &
         5, 6, 4, -1/

      ! K-Ni

      DATA((ncoefs(i, j), i=1, 4), j=19, 28)/ &
         -1, -1, -1, -1, &
         6, 6, 3, -1, &
         -1, -1, -1, -1, &
         -1, -1, -1, -1, &
         -1, -1, -1, -1, &
         6, 2, 4, 3, &
         -1, -1, -1, -1, &
         6, 5, 3, -1, &
         2, 4, 6, 2, &
         4, 4, 2, -1/

      ! Cu-Pt

      DATA((ncoefs(i, j), i=1, 4), j=29, 78)/200*-1/

      ! Au

      DATA(ncoefs(i, 79), i=1, 4)/4, 2, 4, -1/

      ! *****  Edge energies (in eV)  *****

      ! H-Ne

      DATA((edge(i, j), i=1, 3), j=1, 10)/ &
         3*1.d32, &
         3*1.d32, &
         3*1.d32, &
         3*1.d32, &
         3*1.d32, &
         284.0, 2*1.d32, &
         401.0, 2*1.d32, &
         531.7, 2*1.d32, &
         3*1.d32, &
         867.0, 2*1.d32/

      ! Na-Ar

      DATA((edge(i, j), i=1, 3), j=11, 18)/ &
         1071.7, 2*1.d32, &
         49.45, 1303.4, 1.d32, &
         72.78, 1559.9, 1.d32, &
         100.6, 1840.0, 1.d32, &
         3*1.d32, &
         165.0, 2470.5, 1.d32, &
         202.0, 2819.6, 1.d32, &
         245.0, 3202.9, 1.d32/

      ! K-Ni

      DATA((edge(i, j), i=1, 3), j=19, 28)/ &
         3*1.d32, &
         349.31, 4038.1, 1.d32, &
         3*1.d32, &
         3*1.d32, &
         3*1.d32, &
         598.0, 691.0, 5988.8, &
         3*1.d32, &
         707.4, 7111.2, 1.d32, &
         61., 793.8, 7709.5, &
         853.6, 8331.6, 1.d32/

      ! Cu-Pt

      DATA((edge(i, j), i=1, 3), j=29, 78)/150*1.d32/

      ! Au

      DATA(edge(i, 79), i=1, 3)/2220., 2743.9, 1.d32/

      ! Number of edges (-1 => no tabulated data)

      DATA nedges/ &
         0, 0, -1, -1, -1, 1, 1, 1, -1, 1, &
         1, 2, 2, 2, -1, 2, 2, 2, &
         -1, 2, -1, -1, -1, 3, -1, 2, 3, 2, &
         50*-1, 2/

      ! cm*cm/g to barns conversion factor

      DATA f/ &
         1.674, 6.646, 0, 0, 0, 19.94, 23.26, 26.56, 0, 33.50, &
         38.17, 40.35, 44.80, 46.63, 0, 53.24, 58.86, 66.33, 0, 66.54, &
         0, 0, 0, 86.33, 0, 92.72, 97.85, 97.48, 50*0., 327.0/

      photo = 0.

      ! Special case for Helium

      IF (Z .EQ. 2) THEN
         IF (versn .EQ. 2) THEN
            photo = helxsc(500*(keV1 + keV2))*EXP(-55.26204)*f(Z)
         ELSEIF (versn .EQ. 3) THEN
            photo = helyan(500*(keV1 + keV2))*EXP(-55.26204)*f(Z)
         END IF
         RETURN
      END IF

      ! Other elements

      IF (nedges(Z) .EQ. -1) THEN
         status = 1
         photo = 0.
         RETURN
      ELSE
         status = 0
      END IF

      E1 = 1000*keV1

      E2 = 1000*keV2

      ! Find the appropriate range to contain the lower energy.

      i1 = 1
      DO WHILE ((E1 .GE. edge(i1, Z)) .AND. (i1 .LE. nedges(Z)))
         i1 = i1 + 1
      END DO

      ! Find the appropriate range to contain the upper energy.

      i2 = i1
      DO WHILE ((E2 .GE. edge(i2, Z)) .AND. (i2 .LE. nedges(Z)))
         i2 = i2 + 1
      END DO

      ! If these are the same then just sum up the cross-section
      ! for the midpoint of the bin.

      IF (i1 .EQ. i2) THEN

         Elog = log((E1 + E2)/2)

         X = coeff(ncoefs(i1, Z), i1, Z)
         DO j = ncoefs(i1, Z) - 1, 1, -1
            X = X*Elog + coeff(j, i1, Z)
         END DO

      ELSE

         ! First do the lower energy up to the first edge in the bin

         Elog = log((E1 + edge(i1, Z))/2)

         Xt = coeff(ncoefs(i1, Z), i1, Z)
         DO j = ncoefs(i1, Z) - 1, 1, -1
            Xt = Xt*Elog + coeff(j, i1, Z)
         END DO
         X = Xt*(edge(i1, Z) - E1)/(E2 - E1)

         ! Now calculate the last edge in the bin up to the upper energy

         Elog = log((E2 + edge(i2 - 1, Z))/2)

         Xt = coeff(ncoefs(i2, Z), i2, Z)
         DO j = ncoefs(i2, Z) - 1, 1, -1
            Xt = Xt*Elog + coeff(j, i2, Z)
         END DO
         X = X + Xt*(E2 - edge(i2 - 1, Z))/(E2 - E1)

         ! Now add in any bits between edges in the bin

         DO k = i1 + 1, i2 - 1

            Elog = log((edge(k, Z) + edge(k - 1, Z))/2)
            Xt = coeff(ncoefs(k, Z), k, Z)
            DO j = ncoefs(k, Z), 1, -1
               Xt = Xt*Elog + coeff(j, k, Z)
            END DO
            X = X + Xt*(edge(k, Z) - edge(k - 1, Z))/(E2 - E1)

         END DO

      END IF

      ! Do the exponential, put in the E**3 factor, and convert to cm**2.

      Elog = log((E1 + E2)/2)
      photo = EXP(X - 3*Elog - 55.26204)*f(Z)

      RETURN
   END FUNCTION photo

   FUNCTION HELXSC(E)

      REAL HELXSC, E

      !
      !     Real Funcion : HELIUM
      !     Source : Marr, G. V., and West, J. B., Atomic and Nuclear Data Tables,
      !                (1976) 18, 497.
      !             Oza, D. H., (1986), Phys. Rev. A, 33,  824.
      !             Fernley, J. A., Taylor, K. T., and Seaton, M. J., (1987),
      !                J. Phys. B., 20, 6457.
      !
      !     Description :
      !     calculates mass absorption coefficient (mu/rho) in cm2/g for neutral
      !     helium for the given energy in eV.
      !     Cross sections come from experimental data compiled by Marr and
      !     West (Atomic Data and Nuclear Data Tables (1976) 18, 497).
      !     The four strongest autoionization resonances are taken into account;
      !     numbers come from Oza (Phys Rev A (1986), 33, 824), and Fernley et al.
      !     (J. Phys B (1987) 20, 6457).
      !
      !     Deficiencies :
      !     works in the energy range from 30 eV to 10,000 eV

      !     Bugs :
      !     if any are found please report to the authors
      !
      !     History :
      !     this subroutine replaces the previous version of HELIUM which
      !     calculated mass absoprtion coefficients based on Henke's data
      !     (Henke, B. L., et al., (1982), Atomic and Nuclear Data Tables, 27, 1).
      !     This version of HELIUM returns mass  absorption coefficients which
      !     are in better agreement with the best experiments as well as
      !     theoretical models (see Chen, W. F., Cooper, G., and Brion, C. E.,
      !     (1991), Phys. Rev. A, 44, 186).  This fortran-77 version of the
      !     subroutine is based on Pat Jelinsky's program written in C
      !     (obtained from EUVE Archive)
      !
      !     History :
      !     04 jan 93 : original (19775::MBC)
      !
      !     23 feb 93 : comments added and modified to remove VAX
      !                    fortran 77 extensions
      !
      !     21 sep 93 : further remnants of VAX fortran 77 extensions
      !                    have been removed (19775::MBC)
      !
      !     23 sep 93 : bug in the FANO routine has been removed (19775::MBC)
      !
      !     Usage : FUNCTION HELIUM(E)
      !            E = Energy in eV
      !
      !     Common Blocks :
      !           none
      !
      !     Implicit :
      !           none
      !
      !     Functions called by HELIUM
      !           FANO
      !
      !------------------------------------------------------------------------------

      !  Avogadro's number

      REAL AV
      PARAMETER(AV=6.022045E23)

      !  atomic weight of hydrogen

      REAL AW
      PARAMETER(AW=4.0026E0)

      INTEGER IP, IF
      PARAMETER(IP=8, IF=4)

      REAL C1(IP), C2(IP), Q(IF), NU(IF), GAMMA(IF)
      REAL LAMBDA, X, Y, SIGMA, EPS

      INTEGER I

      ! polynomial coefficients for Marr and West data
      DATA C1/-2.953607E1, 7.083061E0, 8.678646E-1, -1.221932E0, &
         4.052997E-2, 1.317109E-1, -3.265795E-2, 2.500933E-3/

      ! polynomial coefficients for Marr and West data )
      DATA C2/-2.465188E1, 4.354679E0, -3.553024E0, 5.573040E0, &
         -5.872938E0, 3.720797E0, -1.226919E0, 1.576657E-1/

      ! parameters Q for resonances (Fernley et al. 1987)
      DATA Q/2.81E0, 2.51E0, 2.45E0, 2.44E0/

      ! parameters NU for resonances (Oza 1986)
      DATA NU/1.610E0, 2.795E0, 3.817E0, 4.824E0/

      ! parameters GAMMA for resonances (Oza 1986)
      DATA GAMMA/2.64061E-3, 6.20116E-4, 2.56061E-4, &
         1.320159E-4/

      ! Calculate wavelength

      LAMBDA = 12398.54E0/E
      X = ALOG10(LAMBDA)

      ! If > 503.97 then no absorption

      IF (LAMBDA .GT. 503.97E0) THEN
         HELXSC = 0.E0
         RETURN

         ! If < 46 then use first polynomial fit

      ELSEIF (LAMBDA .LT. 46.E0) THEN
         Y = 0.E0
         DO I = 1, IP
            Y = Y + C2(I)*(X**(I - 1))
         END DO

         ! Otherwise use second polynomial fit and include autoionization
         ! resonances

      ELSE

         Y = 0.E0
         DO I = 1, IP
            Y = Y + C1(I)*(X**(I - 1))
         END DO

         EPS = 911.2671E0/LAMBDA
         DO I = 1, IF
            X = 2.0*(EPS - 3.0E0 + 1.E0/(NU(I)*NU(I)) - 1.807317)/GAMMA(I)
            Y = Y + ALOG10((X - Q(I))*(X - Q(I))/(1.0E0 + X*X))
         END DO

      END IF

      SIGMA = 10.E0**Y
      HELXSC = SIGMA*AV/AW
   END FUNCTION HELXSC

   FUNCTION HELYAN(E)
      !-----------------------------------------------------------------------------
      !
      !     Real Funcion : HELYAN
      !     Source : Yan et al 1998, ApJ 496, 1044.
      !             Oza, D. H., (1986), Phys. Rev. A, 33,  824.
      !             Fernley, J. A., Taylor, K. T., and Seaton, M. J., (1987),
      !                J. Phys. B., 20, 6457.
      !
      !     Description :
      !     calculates mass absorption coefficient (mu/rho) in cm2/g for neutral
      !     helium for the given energy in eV.
      !     Cross sections come from a theoretical adjustment of experimental data
      !     (from Samson et al 1994a and others) by Yan et al. (ApJ 1997).
      !
      !     The four strongest autoionization resonances are taken into account;
      !     numbers come from Oza (Phys Rev A (1986), 33, 824), and Fernley et al.
      !     (J. Phys B (1987) 20, 6457).
      !
      !     Deficiencies :
      !     works in the energy range from 30 eV to 10,000 eV

      !     Bugs :
      !     if any are found please report to the authors
      !
      !     History :
      !     1991: calculated mass absoprtion coefficients based on Henke's data
      !     (Henke, B. L., et al., (1982), Atomic and Nuclear Data Tables, 27, 1).
      !     1993: modified to include autoionization resonances.  The fortran-77
      !     version of the subroutine was based on Pat Jelinsky's C program.
      !     (obtained from EUVE Archive)
      !     1997: realized Marr and West (1976:
      !     cross sections are a poor match to current best estimates above 80 eV.
      !     Changed continuum portion to cross sections recommended by Yan et al
      !     (ApJ 1977)   Autoionization resonances were
      !     retained as implemented by Jelinsky in the EUVE archive.
      !
      !
      !     Usage : FUNCTION HELYAN(E)
      !            E = Energy in eV
      !
      !     Common Blocks :
      !           none
      !
      !     Implicit :
      !           none
      !
      !     Functions called by HELYAN
      !           FANO
      !
      !------------------------------------------------------------------------------

      !    Type definitions :
      !     IMPLICIT NONE
      !    Global variables :
      !    Structure definitions :
      !    Function declarations :

      !    Local constants :
      INTEGER IP
      !         ( index through loop )
      PARAMETER(IP=6)
      INTEGER IF
      !         ( index through loop )
      PARAMETER(IF=4)
      REAL AV
      !         ( Avogadro's number )
      PARAMETER(AV=6.022045E23)
      REAL AW
      !         ( atomic weight of hydrogen )
      PARAMETER(AW=4.0026E0)
      !    Local variables :
      REAL LAMBDA
      !          ( wavelength in Angstroms)
      REAL X
      REAL Y
      REAL SIGMA
      !          ( cross section in cm2/atom)
      INTEGER I
      !          ( index trough loop)
      !     Import :
      REAL E
      !          ( energy in eV)
      !     Export :
      REAL HELYAN
      !          ( cross section in cm**2/g)
      !    Local data :
      REAL C1(IP)
      REAL EION
      REAL Q(IF)
      REAL NU(IF)
      REAL GAMMA(IF)

      !          ( polynomial coefficients for Yan et al data)
      DATA C1/-4.7416, 14.8200, -30.8678, 37.3584, -23.4585, 5.9133/

      !          ( ionization edge in eV:)
      DATA EION/24.58/

      !          ( parameters Q for resonances (Fernley et al. 1987) )
      DATA Q/2.81E0, 2.51E0, 2.45E0, 2.44E0/

      !          ( parameters NU for resonances (Oza 1986) )
      DATA NU/1.610E0, 2.795E0, 3.817E0, 4.824E0/

      !          ( parameters GAMMA for resonances (Oza 1986) )
      DATA GAMMA/2.64061E-3, 6.20116E-4, 2.56061E-4, &
         1.320159E-4/

      !     Start :

      LAMBDA = 12398.54E0/E

      X = E/EION

      IF (LAMBDA .GT. 503.97E0) THEN
         HELYAN = 0.E0
      ELSE
         Y = 1.E0
         DO 2 I = 1, IP
            Y = Y + C1(I)/(X**(I/2.))
2        CONTINUE

         ! Yan et al. cross section in cm**2/atom:
         SIGMA = 733.E-24/(1.E-3*E)**3.5*Y

         ! Add in autoionization resonances:
         DO 3 I = 1, IF
            SIGMA = SIGMA*FANO(Q(I), NU(I), GAMMA(I), LAMBDA)
3        CONTINUE

         HELYAN = SIGMA*AV/AW

      END IF

      RETURN
   END FUNCTION HELYAN


   FUNCTION FANO(A, B, C, LAMBDA)

      !    Type definitions :
      !     IMPLICIT NONE
      !    Global variables :
      !    Structure definitions :
      !    Function declarations :
      !    Local constants :
      !    Local variables :
      REAL EPS
      !          ( energy in Rydbergs )
      REAL EPSI
      REAL X
      !          ( log_10 of wavelength in Angstroms )
      !     Import :
      REAL A
      !          ( Q coefficient (Fernley et al. 1987) )
      REAL B
      !          ( NU coefficient (Oza 1986) )
      REAL C
      !          ( GAMMA coefficient (Oza 1986) )
      REAL LAMBDA
      !          ( wavelength in Angstroms )
      !     Export :
      REAL FANO
      !    Start :

      EPS = 911.2671E0/LAMBDA
      EPSI = 3.0E0 - 1.E0/(B*B) + 1.807317
      X = 2.0*(EPS - EPSI)/C
      FANO = (X - A)*(X - A)/(1.0E0 + X*X)
   END FUNCTION FANO


   FUNCTION XraySintegrand(zz)

      IMPLICIT NONE

      REAL*8               :: zz, rr
      REAL*8               ::  XraySintegrand

      ! uu is the projected distance from the cluster center on the sky
      ! zz is the distance along the line of sight.
      rr = SQRT(uu*uu + zz*zz)

      ! We compare to rr because we want to match the spherical shape of the cluster

      IF (rr < r_min) THEN
         rr = r_min
         XraySintegrand = Xrayemissfunc1(rr)
      ELSEIF (rr >= r_integration_max) THEN
         XraySintegrand = 0.d0
      ELSE
         XraySintegrand = Xrayemissfunc1(rr)
      END IF

   END FUNCTION XraySintegrand
!=================================================================================================

   FUNCTION Xrayemissfunc1(rr)

      IMPLICIT NONE

      REAL*8               ::rr
      REAL*8               :: Xrayemissfunc1
      REAL*8               :: result

      ! TODO: We're nesting checks in a way that may be unnecessary

      IF (rr < r_min) THEN
         CALL interp1d_even(logX_emiss1D, logr, n, phlog10(r_min), result)

      ELSEIF (rr >= r_integration_max) THEN
         Xrayemissfunc1 = 0.d0
         RETURN
      ELSE
         CALL interp1d_even(logX_emiss1D, logr, n, phlog10(rr), result)
      END IF

      Xrayemissfunc1 = 10.0**result

   END FUNCTION Xrayemissfunc1

!================================================================================================

end MODULE GasModels
