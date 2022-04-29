MODULE params

!========================================================================================================================
!                             Bayes-X.inc
!
! Include file for BAYES-X   M.Olamaie September 2013
!=========================================================================================================================

   INTEGER                                                       :: myID = 0

! Object parameters:

! ??? Are there objects to be fitted?
   INTEGER, PARAMETER                               :: Atoms = 1

! ??? How many?
   INTEGER, PARAMETER                               :: NAtoms = 1
   INTEGER                                            :: tot_atoms
! Available Model Classes:
!     1 = Clusters
!     1 = Gas
   INTEGER, PARAMETER                                 :: ModelClass = 1
   INTEGER, PARAMETER                                 :: Gas = 1
! Available geometry models:
!     1 = Circular symmetry (2 pars, x/arcsec & y/arcsec)
! ??? Which geometry model?
   INTEGER, PARAMETER                                 :: GeoModel = 1
   INTEGER, PARAMETER                                 :: NGeoPars = 2*GeoModel
   REAL*8                                             :: GeoPars(NAtoms, NGeoPars)
! Available gas models:
!     5= DM_GNFW model (6 pars: Mtot(r200), fgas(r200), a, b, c, c500
! ??? Which gas density profile model?
   INTEGER                                              :: GasModel
   INTEGER, PARAMETER                                   :: NGasPars = 7
   REAL*8                                              :: GasPars(NAtoms, NGasPars)
! redshift and angular diameter variables and lookup tables
   INTEGER                                              :: Dn, SCZdn
   REAL*8                                              :: D
   REAL*8                                              :: z(NAtoms)
   REAL*8                                              :: zdmin, zdmax
   REAL*8                                               :: Q(2, 2)
   REAL*8, DIMENSION(:, :), ALLOCATABLE               :: lookM, lookD
   REAL*8, DIMENSION(:, :, :), ALLOCATABLE              :: lookZ

! if znow = T, \rho_{crit} is defined to be the present critical density
! if znow = F, \rho_{crit} is defined to be the critical density at the cluster redshift

   LOGICAL, PARAMETER                                 ::  znow = .FALSE.

! set what mass function to use (relevant only if using joint prior for M200 & z using
!the mass function)
! mass_function = 1, Evrard approximation to Press-Schechter
! mass_function = 2, Jenkins
! mass_function = 3, Tinker

   INTEGER                                              ::  mass_function
   REAL*8, PARAMETER                                          ::  MassLim = 2.d14
   REAL*8, PARAMETER                                          ::  MassMin = 1.d14
   REAL*8, PARAMETER                                          ::  MassMax = 5.d15

!========================================================================================

   INTEGER, PARAMETER                                 ::  nx = 256, ny = 256
   INTEGER                                              ::  n
   INTEGER, PARAMETER                                  ::      aux_dim = 136
   REAL*8                                               ::  aux(NAtoms, aux_dim)
   REAL*8, DIMENSION(:), ALLOCATABLE               ::  r, logr
   REAL*8                                               ::  uu, loguu
   REAL*8, PARAMETER                              ::rmin = 0.01, rmax = 10.0, rlimit = 10.0

   REAL*8                                                ::  rhocritz

! cluster data and working arrays for DM_GNFW model

   REAL*8               ::  a_GNFW, b_GNFW, c_GNFW, c500_GNFW
   REAL*8               ::  rs_DM, rhos_DM, Pei_GNFW, Pei_GNFW_keV, rp_GNFW
   REAL*8               ::  alpha_Einasto, r_2_DM, rho_2_DM
   REAL*8               ::  Gamma_coeff1, Gamma_coeff2, mass_coeff_Einasto
   REAL*8               ::  r200_DM, c200_DM
   REAL*8               ::  Mg200_DM, MT200_DM, fg200_DM
   REAL*8               ::  Tg200_DM, n_e200, Ke200, Pe200
   REAL*8               ::  n_H200, ne_nH200, Rhogas200
   REAL*8               ::  r500_DM, c500_DM
   REAL*8               ::  Mg500_DM, MT500_DM, fg500_DM
   REAL*8               ::  Tg500_DM, n_e500, Ke500, Pe500
   REAL*8               ::          n_H500, ne_nH500, Rhogas500
   REAL*8               ::   r2500_DM, c2500_DM
   REAL*8               ::   Mg2500_DM, MT2500_DM, fg2500_DM
   REAL*8               ::   Tg2500_DM, n_e2500, Ke2500, Pe2500
   REAL*8               ::          n_H2500, ne_nH2500, Rhogas2500
   REAL*8               ::  rx_incre

   REAL*8, DIMENSION(:), ALLOCATABLE  :: rgx, rhogasx, n_Hx, ne_nHx
   REAL*8, DIMENSION(:), ALLOCATABLE  :: n_ex, Tgx, Kex, Pex, M_DMx, Mg_DMx, fg_DMx
!==========================================================================================================================================
! XRAY variables , parameters and working arrays

! Info for MEKAL model:

   INTEGER NUMION
   PARAMETER(NUMION=207)

   INTEGER NREC, NAMAX
   PARAMETER(NREC=20, NAMAX=NREC*256)

   INTEGER NSTO
   PARAMETER(NSTO=5500)

   INTEGER L_CAFE
   PARAMETER(L_CAFE=400)

   INTEGER ne, ncafe

!       ne        from Fillin
!       ncafe     from Filcaf

   INTEGER NUU, NG
   PARAMETER(NUU=25, NG=25)

!       from Filrec...
   REAL arecs(13, NUMION)
   REAL*8 arec(13, NUMION)

!       from Filion...
   REAL aions(NAMAX)
   REAL*8 aion(NAMAX)

!       from Fillin

   INTEGER*2   :: nrl(NSTO)
   INTEGER*2   :: lp(9, NSTO)
   CHARACTER*8 :: trans(NSTO)
   REAL        :: rpsing(14, NSTO)
   REAL*8     :: rp(14, NSTO)

!       from Filcaf

   INTEGER*4   ::  idnr(L_CAFE)
   REAL        ::  cdrcafes(L_CAFE)
   REAL*8     ::  cdrcafe(L_CAFE)

!       from Filspe
   REAL         :: psing(6, NUMION)
   REAL*8      :: p(6, NUMION)

!       from Filkar
   REAL         :: gasing(NG, NUU)
   REAL*8      :: ga(NG, NUU)

! End of info for MEKAL model

   LOGICAL         ::  qanyf
   REAL*8          ::  Rhogas0, logRhogas0, KB_Tx, n_e0, n_H0, mu_ex, mu_mx
   REAL*8          ::  Rhogas0_SI, el_den, sum_xenuctot, sum_xeZtot
   REAL*8          ::  ne_nH0
   INTEGER         ::  indexy, init, lun, ierr, nl, nlx
   CHARACTER(LEN=256)   ::  filion = '', filrec = '', filkar = '', filcon = '', fillin = ''
   CHARACTER(LEN=256)   :: filcaf = ''
   CHARACTER(LEN=256)   ::  filARF = '', filRMF = '', filevent = '', filBG = ''
   CHARACTER(LEN=72)    ::  comment
   INTEGER, PARAMETER  ::  NOEL = 15
   INTEGER, PARAMETER  ::  NL_MAX = 5500

   CHARACTER(LEN=3), PARAMETER   ::  Elements(1:15) = &
                                    (/'H  ', 'He ', 'C  ', 'N  ', 'O  ', 'Ne ', 'Na ', 'Mg ', 'Al ', &
                                      'Si ', 'S  ', 'Ar ', 'Ca ', 'Fe ', 'Ni '/)
   REAL*8              ::  abund, nuc_tot, Z_tot, xe, xzin, elx, flx
   DIMENSION abund(NOEL), nuc_tot(NOEL), Z_tot(NOEL), xe(NOEL), xzin(NUMION), elx(NL_MAX)
   DIMENSION flx(NL_MAX)
   DATA abund/12.00, 10.99, 8.56, 8.05, 8.93, 8.09, 6.33, 7.58, 6.47, 7.55, 7.21, 6.56, 6.36, 7.67, 6.25/
   DATA nuc_tot/1.0, 4.0, 12.0, 14.0, 16.0, 20.0, 23.0, 24.0, 27.0, 28.0, 32.0, 40.0, 40.0, 56.0, 59.0/
   DATA Z_tot/1.0, 2.0, 6.0, 7.0, 8.0, 10.0, 11.0, 12.0, 13.0, 14.0, 16.0, 18.0, 20.0, 26.0, 28.0/

   DATA init/0/

   CHARACTER(LEN=3)                 :: XrayTelescope
   INTEGER                          :: xraynx, xrayny
   INTEGER                          :: xraynxcentre, xraynycentre
   REAL*8                          :: Aeffave
   REAL*8                          :: xraycell
   INTEGER                          :: xrayxpix, xrayypix
   REAL*8                          :: xrayx0, xrayy0
   REAL*8                          :: xrayxx, xrayyy
   REAL*8                          :: xrayr
   REAL*8                          :: xraydx(2)
   REAL*8                          :: xraytrans(6)
   INTEGER                          :: xrayNbin
   INTEGER                          ::  xrayNch
   INTEGER                          :: LENx
   INTEGER                             :: xLENi
   REAL*8                          :: xrayEmin
   REAL*8                          :: xrayEmax
   REAL*8                          :: xrayDeltaE
   REAL*8                          :: sexpotime
   REAL*8                          :: bexpotime
   REAL*8                          :: xrayde
   REAL*8                          :: xfluxsec1, xfluxsec2, xfluxsec3
   REAL*8                           :: xfluxsec4, xfluxsec5, xfluxsec6
   REAL*8, DIMENSION(:), ALLOCATABLE     ::n_H, ne_nH, n_e, T, Rhogas
   REAL*8, DIMENSION(:, :), ALLOCATABLE     :: X_emiss2D, X_S2D, predX_S2D
   REAL*8, DIMENSION(:), ALLOCATABLE              :: logX_emiss1D, X_S1D
   REAL*8, DIMENSION(:), ALLOCATABLE     :: xrayE1, xrayE2
   REAL*8, DIMENSION(:), ALLOCATABLE     :: xrayFluxCoeff, xrayFlux, xraypredFlux
   REAL*8, DIMENSION(:), ALLOCATABLE      :: Arf
   REAL*8, DIMENSION(:, :), ALLOCATABLE     :: Rmf, TRM
   REAL*8, DIMENSION(:), ALLOCATABLE       :: xrayBG
   REAL*8, DIMENSION(:, :, :), ALLOCATABLE   :: xrayCmap
   REAL*8, DIMENSION(:), ALLOCATABLE :: xrayCpred
   INTEGER, DIMENSION(:), ALLOCATABLE :: xrayCobs, xrayBG_obs
   REAL*8                    :: XRAYLhood0
!Hydrogen Column Density
   REAL                    ::  N_H_col
!Background model
   REAL*8                 ::  xrayBG_predmax

!Photoelectric absorption
   REAL, DIMENSION(:), ALLOCATABLE       ::  photar
!=========================================================================================================================

! Priors for each parameter:
! Available Priors:
!     0 = Delta function prior (parameter fixed) - (fixedpt,*)
!     1 = Uniform prior - (xmin,xmax)
!     2 = Uniform prior in log_10 - (xmin,xmax)
!     3 = Gaussian prior - (xmean,xstdev)
!     4 = LogNormal prior - (xmean,xwidth)
!     5 = Sinusoidal prior, uniform in cos - (thetamin,thetamax)
!     6 = Cauchy prior - (xmean,xwidth)
!     7 = Prior for spectral indices from Waldram 07
!     12 = Prior for spectral indices from Waldram 07 multipled by a Gaussian
!     8 = joint prior for M200 & z using the mass function (with mass_function = 1, 2 & 3 for Evrard, Jenkins & Tinker mass functions respectively)
!     9 = joint prior on x & y with (x,y) lying uniformly in a triangle defined by the 3 vertices
!     10 = Exponential prior = lambda * exp(-lambda * x) (min, max, lambda)
!     11 = Power law prior = x^{-alpha} (min, max, alpha)
!
!  XXX_Prior(i,j,k) -
!       i = atom number
!       j = parameter number
!       k = 1,2 - prior definition

   INTEGER          :: Geo_PriorType(NAtoms, NGeoPars)
   REAL*8           :: Geo_Prior(NAtoms, NGeoPars, 3)
   INTEGER          :: Gas_PriorType(NAtoms, NGasPars)
   REAL*8           :: Gas_Prior(NAtoms, NGasPars, 3)
   INTEGER          :: z_PriorType(NAtoms)
   REAL*8           :: z_Prior(NAtoms, 3)

!========================================================================================
!Number of parameters per object

   INTEGER, PARAMETER           :: NPars = NGeoPars + NGasPars*Gas + 1
   INTEGER                       :: NPars2

!Total no. of dimensions of parameter space:
   INTEGER, PARAMETER           :: NDim = NAtoms*NPars
   INTEGER                       :: tot_dim, edim

!==========================================================================================================
! ??? Debugging (may not be very helpful...):
   LOGICAL, PARAMETER          :: verbose = .FALSE.

!whether to do multimodal sampling

   LOGICAL                  :: n_mmodal = .FALSE.
   LOGICAL                         :: n_ceff = .FALSE.
   INTEGER                         :: n_nlive = 1000
   REAL*8                   :: n_efr = 0.8
   INTEGER                         :: n_rseed = -1
   INTEGER                         :: n_totPar
   LOGICAL                         :: n_IS = .false.
   INTEGER                         :: n_updint = 100
   REAL*8                   :: n_tol = 0.5
   CHARACTER(LEN=1000)      :: n_root = ''
   INTEGER, DIMENSION(:), ALLOCATABLE        :: n_pWrap
   LOGICAL                                        :: n_fb = .TRUE.
   INTEGER                                        :: n_maxModes = 20

   REAL*8, PARAMETER                          :: B_A_D = -1.0d64

!=======================================================================

END MODULE params
