MODULE MassModels
	USE params
	USE constants
	USE matrix_utils
	USE utilities
	USE cosmology

CONTAINS

!=======================================================================


      SUBROUTINE EllGeometry(i)

            IMPLICIT NONE

               INTEGER                         ::i	
               REAL*8                          :: theta,f,RR(2,2),SS(2,2),temp1(2,2),temp2(2,2)

!-----------------------------------------------------------------------

! Construct the transformation matrix Q for elliptical geometry fo
! the ith atom:

!   theta=orientation angle of major axis anticlockwise from x axis
!   f=axis ratio b/a < 1

                   theta=GeoPars(i,3)*Pi/180.d0
                   f=GeoPars(i,4)

                   RR(1,1)=cos(theta)	  
                   RR(1,2)=sin(theta)	  
                   RR(2,1)=-sin(theta)	  
                   RR(2,2)=cos(theta)	  
                   SS(1,1)=sqrt(f)	  
                   SS(1,2)=0.0	  
                   SS(2,1)=0.0	  
                   SS(2,2)=1.0/sqrt(f)	  
  
                   CALL MatrixProduct(2,2,SS,2,2,RR,temp1)
                   CALL Transpose(2,2,temp1,temp2)
                   CALL MatrixProduct(2,2,temp2,2,2,temp1,Q)

      END SUBROUTINE EllGeometry

!-----------------------------------------------------------------------
      SUBROUTINE makeDlookup(zdmin,zdmax)

            IMPLICIT NONE

               INTEGER                         :: i,j,k
               REAL*8                          :: zdmin,zdmax

!iterate over lens redshift
                   DO i=1,Dn
                      IF(i==1) THEN
                         lookD(i,1)=zdmin
                      ELSEIF(i==Dn) THEN
                         lookD(i,1)=zdmax
                      ELSE
                         lookD(i,1)=10.**(log10(zdmin)+(log10(zdmax)-log10(zdmin))*(i-1.)/(Dn-1.))
                      ENDIF
                         lookD(i,2)=angdist(lookD(i,1))
                   ENDDO

      END SUBROUTINE makeDlookup

!-----------------------------------------------------------------------
END MODULE MassModels
