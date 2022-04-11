MODULE cosmology

	USE constants
        USE utilities
      
CONTAINS

!=======================================================================

      FUNCTION Hubble(zloc)

      	IMPLICIT NONE
      	
      	REAL*8                        ::  Hubble,zloc,Q

      	IF (wa /= 0.) THEN
        	Q = (1.+zloc)**(3.*(w0+wa))*exp(-3.*wa*zloc/(1.+zloc))
      	ELSE
        	Q = (1.+zloc)**(3.*w0)
      	ENDIF
	
! 	General expression:
      
      	Hubble = (1.+zloc)*(1.+zloc)*(1.+zloc)*(Om + OL*Q) + Ok*(1.+zloc)*(1.+zloc)
      	Hubble = h*100.*sqrt(Hubble)

	END FUNCTION Hubble
     
!=======================================================================
! Function returns rhocrit(zz)

      FUNCTION rhocritofz(zz)

      	IMPLICIT NONE
      
      	REAL*8                        ::  rhocritofz,zz
      
      	IF (zz == 0.) THEN
        	rhocritofz = rhocrit
      	ELSE
        	rhocritofz = rhocrit*(Hubble(zz)/(h*100.))**2.
	ENDIF
            
      END FUNCTION rhocritofz

!=======================================================================
! Function returns rhobar(zz)

      FUNCTION rhobarofz(zz)

      	IMPLICIT NONE
      
      	REAL*8                        ::  rhobarofz,zz
      
      	rhobarofz = Omofz(Om,zz)*rhocritofz(zz)
            
      END FUNCTION rhobarofz

!=======================================================================
! Function returns Omega_m(zz)

      FUNCTION Omofz(Omega_m,zz)

      	IMPLICIT NONE
      
      	REAL*8                        ::  Omofz,Omega_m,zz
      
      	IF (zz == 0.) THEN
        	Omofz = Omega_m
      	ELSE  
        	Omofz = Omega_m*(1.+zz)**3./(1.-Omega_m+Omega_m*(1.+zz)**3.)
	ENDIF
            
      END FUNCTION Omofz

!=======================================================================
! Function the (angular diameter distance) * h as a function of redshift z

      FUNCTION calcAngDiamDis(z1,z2)

      	IMPLICIT NONE
	
      	REAL*8                                     ::  calcAngDiamDis,z1,z2,S,lim1,lim2
      	REAL*8 , PARAMETER                         ::  eps=1d-3
      
      	lim1=1./(1.+z2)
      	lim2=1./(1.+z1)
      	CALL qtrap(AngDiamInt,lim1,lim2,eps,S)
      	calcAngDiamDis=clight*S/(h*100000.*(1.+z2))

	END FUNCTION calcAngDiamDis
	
!=======================================================================
      
      FUNCTION AngDiamInt(r)
      
      	IMPLICIT NONE
      
      	REAL*8                        ::  r,AngDiamInt
      
      	AngDiamInt=1./sqrt(Om*r+OL*r**(1.-3.*w0))
      
      END FUNCTION AngDiamInt
	
!=======================================================================
! Comoving distance in h^-1 Mpc:

	FUNCTION rcomoving(zz)
      
      	IMPLICIT NONE
      
      	REAL*8                                   ::  rcomoving,zz
      	REAL*8 , PARAMETER                       ::  eps=1d-4
      
      	CALL qtrap(conH,0.d0,zz,eps,rcomoving)
            
      END FUNCTION rcomoving

!=======================================================================

	FUNCTION conH(zz)
      
      	IMPLICIT NONE
      
      	REAL*8                                   :: conH,zz
      
      	IF (zz == 0.0) THEN
        	conH = 2.998d5/100.0
      	ELSE
        	conH = 2.998d5/Hubble(zz)
		ENDIF
            
      END FUNCTION conH

!=======================================================================
! Return lookback time in Gyr:
	FUNCTION lookbacktime(zz)
      
     	 	IMPLICIT NONE
      
      	REAL*8                        ::  lookbacktime,zz
      	REAL*8 , PARAMETER            ::  eps=1d-3
      	
      	CALL qtrap(timeintegrand,0.d0,zz,eps,lookbacktime)
      
      END FUNCTION lookbacktime

!=======================================================================

	FUNCTION timeintegrand(zz)
      
      	IMPLICIT NONE
      
      	REAL*8                        ::  timeintegrand,zz
      
      	timeintegrand = 100./(0.1021*Hubble(zz)*(1.+zz))
      
      END FUNCTION timeintegrand

!=======================================================================

	FUNCTION Scurvature(r)
      
      	IMPLICIT NONE
      
      	REAL*8                         ::  Scurvature,r
      	INTEGER                        ::  k
      
		IF (Ok == 0.) THEN
	  		k = 0
		ELSEIF (Ok > 0.) THEN
	  		k = -1
		ELSEIF (Ok < 0.) THEN
	  		k = 1
		ENDIF 
      
		IF (k == -1) THEN
	  		Scurvature = sinh(r)
		ELSEIF (k == 0) THEN
	  		Scurvature = r
		ELSEIF (k == 1) THEN
	  		Scurvature = sin(r)
		ENDIF  
      
      END FUNCTION Scurvature

!=======================================================================

	FUNCTION angdist(zz)
      
      	IMPLICIT NONE
      
      	REAL*8                         :: angdist,zz
      	REAL*8                         ::  r,R0
      
	      IF (Ok == 0.) THEN
      	  	R0 = conH(0.d0)
      	ELSE
        	R0 = conH(0.d0)/sqrt(abs(Ok))
      	ENDIF
      	r = rcomoving(zz)
      	angdist = R0*Scurvature(r/R0)/(1.+zz)
      
      END FUNCTION angdist

!=======================================================================

	FUNCTION angdistdiff(z1,z2)
      
      	IMPLICIT NONE
      
	REAL*8                        ::  angdistdiff,z1,z2
        REAL*8                        ::  r1,r2,dr,R0
      
	      IF (Ok == 0.0) THEN
        		R0 = conH(0.d0)
      	ELSE
        		R0 = conH(0.d0)/sqrt(abs(Ok))
      	ENDIF
      	r1 = rcomoving(z1)
      	r2 = rcomoving(z2)
      	angdistdiff = R0*Scurvature((r2-r1)/R0)/(1.+z2)
      
      END FUNCTION angdistdiff

!=======================================================================

	FUNCTION criticaldensity(z1,z2)
      
      	IMPLICIT NONE
            
      	REAL*8                        ::  criticaldensity,z1,z2
      
      	criticaldensity = angdist(z2)/(angdist(z1)*angdistdiff(z1,z2))
      	criticaldensity = criticaldensity*1.6620d6
      
      END FUNCTION criticaldensity

!=======================================================================

	FUNCTION lumdist(zz)
      
      	IMPLICIT NONE
      
      	REAL*8                        ::  lumdist,zz
      	REAL*8                         :: R0
      
      	IF (Ok == 0.) THEN
        		R0 = conH(0.d0)
      	ELSE
        		R0 = conH(0.d0)/sqrt(abs(Ok))
      	ENDIF
      	lumdist = R0*Scurvature(rcomoving(zz)/R0)*(1.+zz)
      
      END FUNCTION lumdist

!=======================================================================
! 	Differential comoving volume for unit solid angle:

	FUNCTION dVdzcomoving(zz)
      
      	IMPLICIT NONE
      
      	REAL*8                        ::  dVdzcomoving,zz
      	REAL*8                        ::  S,r,R0
      
      	IF (Ok == 0.0) THEN
        		R0 = conH(0.d0)
      	ELSE
        		R0 = conH(0.d0)/sqrt(abs(Ok))
      	ENDIF
      	r = rcomoving(zz)
      	S = R0*Scurvature(r/R0)
      	dVdzcomoving = conH(zz)*S*S
      
      END FUNCTION dVdzcomoving

!=======================================================================
! 	Comoving volume between 2 redshifts for unit solid angle:

	subroutine Vcomoving(z1,z2,Volume)
      
      	IMPLICIT NONE
      
      	REAL*8                        ::  z1,z2,Volume
      	REAL*8 , PARAMETER            ::  eps=1d-3
      	
      	CALL qtrap(dVdzcomoving,z1,z2,eps,Volume)
      
      END subroutine Vcomoving

!=======================================================================

END MODULE cosmology
