MODULE CheckPars1
	USE params
	USE constants
	USE utilities
	USE matrix_utils
	USE MassModels
	USE priors1
	USE massfunction
        USE GasModels 

	    IMPLICIT NONE
	
	       LOGICAL crapAlpha
	       LOGICAL :: noCluster=.FALSE.

CONTAINS
!=======================================================================

	SUBROUTINE CheckPars(Cube,flag)

	    IMPLICIT NONE
	
	        REAL*8                         ::  Cube(Ndim)
	        INTEGER flag,j

! Rescale parameters:

	            DO j = 1, Ndim
	                 crapAlpha = .FALSE.
	  		 CALL Rescale(Cube(j),j)
		         IF( crapAlpha ) THEN
	                     flag = 1
	                     RETURN
	                 ENDIF
	            ENDDO
  	            CALL CheckzPars(flag)
      
END SUBROUTINE CheckPars
	
!=======================================================================


	SUBROUTINE CheckzPars(flag)
      
	        INTEGER                         :: flag, k

		flag = 0      
	            DO k = 1, Natoms
	                IF( z_PriorType(k) > 8 .OR. z_PriorType(k) == 5 .OR. z_PriorType(k) == 6 .OR. z_PriorType(k) == 7 ) THEN
	                    if( myID == 0 ) WRITE(*,*)"ERROR: wrong prior class used for redshift"
	                    flag = 1
	                ENDIF
	            ENDDO
      
	END SUBROUTINE CheckzPars

!=======================================================================	

	SUBROUTINE Rescale(value,ival)
	
	    IMPLICIT NONE

	        INTEGER                         ::  ival,i,ii,j,k,ih,jh,tot_pars
	        REAL*8                          ::  value,dummy
	        REAL*8                          ::  a(2), b(2), c(2), p(2)

	            tot_pars=NPars
		
!     	Rescale the parameters of the atom, labelled j=1..NAtoms; i runs from 1 to NPars:

	            j=(ival-1)/NPars+1
	            i=ival-(j-1)*NPars
	
	            IF(i <= NGeoPars) THEN	
	               ii = i
		       IF( Geo_PriorType(j,ii) == 9 .AND. ii <= 2 ) THEN
	                   GeoPars(j,ii)=Prior(1,value,Geo_Tri_Prior(j,ii,1),Geo_Tri_Prior(j,ii,2),Geo_Tri_Prior(j,ii,3))
                            !cluster position in triangle constraint
	                    IF( ii == 2 ) THEN
	                        a(1:2) = Geo_Prior(j,1:2,1)
	                        b(1:2) = Geo_Prior(j,1:2,2)
	                        c(1:2) = Geo_Prior(j,1:2,3)
	                        p(1:2) = GeoPars(j,1:2)
	                        IF( .not.inTriangle(a, b, c, p) ) crapAlpha = .TRUE.
	                    ENDIF
	               ELSE
	                    GeoPars(j,ii)=Prior(Geo_PriorType(j,ii),value,Geo_Prior(j,ii,1),Geo_Prior(j,ii,2),Geo_Prior(j,ii,3))
	               ENDIF
	            ELSEIF( i <= NGeoPars + 1 ) THEN
	                    IF(z_PriorType(j)==8) THEN
                               !mass function prior THEN call the lookup table
	                       z(j) = value
                               !CALL lookUpZ( value, MassPars(j,2), z(j) )
	                    ELSE
	                        z(j) = Prior(z_PriorType(j), value, z_Prior(j,1), z_Prior(j,2),z_Prior(j,3))
	                        IF( z(j) < 0.01d0 .OR. z(j) > zdmax ) crapAlpha = .TRUE.
	                    ENDIF  		
	            ELSEIF( i <= ( NGeoPars + 1 + NGasPars * Gas ) ) THEN
	                    ii = i - NGeoPars - 1 
	                     IF(  ii==1.AND. Gas_PriorType(j,ii)==8) THEN
                                  !mass function prior for M200
                                  !call the lookup table
	                          CALL lookUpM(value, GasPars(j,ii))
	                          IF(z_PriorType(j)==8) THEN
	                               CALL lookUpZ(z(j), GasPars(j,ii), z(j))
	                               IF( z(j) < 0.01d0 .OR. z(j) > zdmax ) crapAlpha = .TRUE.
	                          ENDIF
	                     ELSE	 
	                          GasPars(j,ii)=Prior(Gas_PriorType(j,ii),value,Gas_Prior(j,ii,1),Gas_Prior(j,ii,2),Gas_Prior(j,ii,3))
	                     ENDIF    	
	            ENDIF

	
	END SUBROUTINE Rescale
	
!=======================================================================

	SUBROUTINE Rescale_nest(value,ival1,cont)
	
	    IMPLICIT NONE
	
	        INTEGER                         ::  ival1,ival,i,ii,j,ih,jh,tot_pars,k,cont
	        REAL*8                          :: value,dummy

	            ival=ival1
	            tot_pars=NPars

!     	Rescale the parameters of the atom, labelled j=1..NAtoms; i runs from 1 to NPars:

	            j=(ival-1)/NPars+1
	            i=ival-(j-1)*NPars
	            IF(i <= NGeoPars) THEN
	                ii = i
	                value = GeoPars(j,ii)
	            ELSEIF( i <= NGeoPars + 1 ) THEN
	                value = z(j)                      	
	            ELSEIF( i <= ( NGeoPars + 1 +  NGasPars * Gas ) ) THEN
	                ii = i - NGeoPars - 1
	            	value = GasPars(j,ii)
	            ENDIF


	END SUBROUTINE Rescale_nest

!=======================================================================

       SUBROUTINE PClass(ival,Ptype)
	
	    IMPLICIT NONE
	        INTEGER                         ::  ival,i,ii,j,ih,jh,tot_pars	
	        INTEGER                         ::  Ptype		
        	

	            tot_pars=NPars
	   	    j = ( ival -1 ) / tot_pars + 1
	   	    i = ival - ( j - 1 ) * tot_pars
         	
            	    IF( i <= NGeoPars )THEN	
         		   ii = i	 
			   Ptype = Geo_PriorType(j,ii)
            	    ELSEIF( i <= NGeoPars + 1 )THEN
            		   Ptype = z_PriorType(j)
         	
         	   ELSEIF(  i <= ( NGeoPars + 1 +  NGasPars * Gas ) ) THEN
			        ii = i - NGeoPars - 1 	 
			        Ptype = Gas_PriorType(j,ii)
      		
	   	    ENDIF
      	

       END SUBROUTINE PClass
	
!=======================================================================


END MODULE CheckPars1
