MODULE constants

	USE params
	
	IMPLICIT NONE
	
! Global variables
	REAL*8                                   ::  Geo_Tri_Prior(NAtoms, 2, 3)
	INTEGER, ALLOCATABLE                     ::  zsrcID(:)
! Physical constants:
	REAL*8,PARAMETER                         :: TCMB=2.726
	REAL*8,PARAMETER                         :: hplanck=6.6260755d-34
	REAL*8 ,PARAMETER                        :: kboltzmann=1.380658d-23
	REAL*8 ,PARAMETER                        :: clight=299792458.
	REAL*8 ,PARAMETER                        :: planckfactor=hplanck/(kboltzmann*TCMB)
        REAL*8 ,PARAMETER                        :: sigma_T=6.6524d-29	
        REAL*8 ,PARAMETER                        :: mec2=511.d0	
        REAL*8,PARAMETER                         :: mu_e=1.14*1.67262158d-27	
        REAL*8,PARAMETER                         :: mu_m=0.6*1.67262158d-27	
        REAL*8 ,PARAMETER                        :: m_sun=1.98892d+30	
! Numerical constants:
      REAL*8 ,PARAMETER                          :: Pi= 3.1415926535d0    
      REAL*8 ,PARAMETER                          :: SqrtPi= 1.772453851d0   
      REAL*8 ,PARAMETER                          :: TwoPi = 6.283185307d0      
      REAL*8 ,PARAMETER                          :: SqrtTwoPi= 2.506628275d0     
! Conversion factors
      REAL*8 ,PARAMETER                          :: sec2rad=4.8481368d-6
      REAL*8 ,PARAMETER                          :: rad2sec=1.0/sec2rad      
      REAL*8 ,PARAMETER                          :: min2rad=sec2rad*60.0
      REAL*8 ,PARAMETER                          :: rad2min=1.0/min2rad      
      REAL*8 ,PARAMETER                          :: deg2rad=sec2rad*3600.0
      REAL*8 ,PARAMETER                          :: rad2deg=1.0/deg2rad     
      REAL*8 ,PARAMETER                          :: sec2min =1./60.         
      REAL*8 ,PARAMETER                          :: min2sec=60.d0      
      REAL*8 ,PARAMETER                          :: Mpc2m=3.08568025d22 !mega parsec to metres      
      REAL*8,PARAMETER                           :: m2Mpc=3.24077649d-23 !metres to mega parsec           
      REAL*8 ,PARAMETER                          :: m2cm=1.0d2    
      REAL*8 ,PARAMETER                          :: J2keV =6.24251d+15          	
      REAL*8 ,PARAMETER                          :: SigmaSq2G = 1.154d8    
      REAL*8 ,PARAMETER                          :: Gmu=3.12d-14
      REAL*8 ,PARAMETER                          :: c_mpc=9.7156d-15 !speed of light in Mpc/s
      REAL*8 ,PARAMETER                          :: G=4.518d-48 !Gravitational constant in Mpc^3 M_sun^-1 s^-2   
      REAL*8 ,PARAMETER                          :: m_p=8.40969762d-58 !proton mass in solar masses     
      REAL*8 ,PARAMETER                          :: m_p_kg=1.672622d-27 !proton mass in kg            
      REAL*8 ,PARAMETER                          :: sigma_T_Mpc2=sigma_T/(Mpc2m * Mpc2m)         
     REAL*8 ,PARAMETER                           :: mec2_keV=511.11

! Cosmology

	REAL*8 ,PARAMETER                       :: h=0.7
	REAL*8 ,PARAMETER                       :: Om =0.30
	REAL*8 ,PARAMETER                       :: OL =0.70	
	REAL*8 ,PARAMETER                       :: Ok=1.-Om-OL
	REAL*8 ,PARAMETER                       :: Ob=0.041
	REAL*8 ,PARAMETER                       :: w0=-1.
	REAL*8 ,PARAMETER                       :: wa=0.
	REAL*8 ,PARAMETER                       :: sigma8=0.8
	REAL*8 ,PARAMETER                       :: rhocrit=3.*(h*100./3.0857e19)**2./(8.*Pi*G)

   

END MODULE constants
