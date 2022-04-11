MODULE  ReadWrite
      USE params
      USE utilities
      USE constants
      USE MassModels

CONTAINS

       SUBROUTINE  ReadInData
	
            IMPLICIT NONE
	
	        INTEGER                         ::  i,j,m,iend
	        INTEGER                         ::  ii,jj,mm	
	        REAL*8                          ::  ra,dec,ra2,dec2,rao,deco
	        CHARACTER(LEN=100)    :: string

         ! Reading X ray data
      
               ! First read data for MEKAL Model

                    IF ( init.EQ.0 ) THEN
                          init = 1
	                 CALL ftopen(lun,filion,0,indexy,ierr)
                         CALL FTMRHD(lun,1,indexy,ierr)
                         DO  i = 1 , NREC
                              CALL FTGCVE(lun,1,i,1,256,0.,AIOns((i-1)*256+1),qanyf,ierr)
                         ENDDO
                         CALL FTCLOS(lun,ierr)
                         IF ( ierr .NE. 0 ) THEN
                           if( myID == 0 ) WRITE(*,*) 'ERROR: Failed to get data from mekal1.dat'
                           STOP
                         ENDIF
                         DO i=1,NAMAX
	                      aion(i)=   AIOns(i)
                         ENDDO
                         CALL ftopen(lun,filrec,0,indexy,ierr)
                         CALL FTMRHD(lun,1,indexy,ierr)
                         DO  i = 1 , NUMION
                               CALL FTGCVE(lun,1,i,1,13,0.,AREcs(1,i),qanyf,ierr)
                         ENDDO
                         CALL FTCLOS(lun,ierr)
                         IF ( ierr .NE. 0 ) THEN
                             if( myID == 0 ) WRITE(*,*) 'Failed to get data from mekal2.dat'
                             STOP
                         ENDIF
                         DO i=1,13
	                      DO j=1,NUMION
	                            arec(i,j)=AREcs(i,j)
	                      ENDDO
	                 ENDDO
                         CALL ftopen(lun,filkar,0,indexy,ierr)
	                 CALL FTMRHD(lun,1,indexy,ierr)
                         DO  i = 1 , NUU
                               CALL FTGCVE(lun,1,i,1,NG,0.,GAsing(1,i),qanyf,ierr)
                         ENDDO
                         CALL FTCLOS(lun,ierr)
                         IF ( ierr .NE. 0 ) THEN
                             if( myID == 0 ) WRITE(*,*) 'Failed to get data from mekal3.dat' 
                             STOP
                         ENDIF
                         DO i=1,NG
	                       DO j=1,NUU
	                            ga(i,j)=gAsing(i,j)
	                        ENDDO
	                 ENDDO	 
                         CALL ftopen(lun,filcon,0,indexy,ierr)
                         CALL FTMRHD(lun,1,indexy,ierr)
                         DO  i = 1 , NUMION
                               CALL FTGCVE(lun,1,i,1,6,0.,Psing(1,i),qanyf,ierr)
                          ENDDO
                         CALL FTCLOS(lun,ierr)
                         IF ( ierr .NE. 0 ) THEN
                              if( myID == 0 ) WRITE(*,*) 'Failed to get data from mekal4.dat' 
                              STOP
                         ENDIF
                         DO i=1,6
	                      DO j=1,NUMION
	                           p(i,j)=Psing(i,j)
	                      ENDDO
	                 ENDDO	 	 
                         CALL ftopen(lun,fillin,0,indexy,ierr)
                         CALL FTMRHD(lun,1,indexy,ierr)
                         CALL FTGKYJ(lun,'NAXIS2',NE,comment,ierr)
                         DO  i = 1 , NE
                               CALL FTGCVI(lun,1,i,1,1,0,NRL(i),qanyf,ierr)
                               CALL FTGCVI(lun,2,i,1,9,0,LP(1,i),qanyf,ierr)
                               CALL FTGCVS(lun,3,i,1,1,' ',TRAns(i),qanyf,ierr)
                               CALL FTGCVE(lun,4,i,1,14,0.,RPsing(1,i),qanyf,ierr)
                         ENDDO
                         CALL FTCLOS(lun,ierr)
                         IF ( ierr .NE. 0 ) THEN
                             if( myID == 0 ) WRITE(*,*) 'Failed to get data from mekal5.dat' 
                             STOP
                         ENDIF
                         DO i=1,14
	                      DO j=1,NSTO
	                           rp(i,j)=RPsing(i,j)
	                      ENDDO
	                 ENDDO	 		 
                         CALL ftopen(lun,filcaf,0,indexy,ierr)
                         CALL FTMRHD(lun,1,indexy,ierr)
                         CALL FTGKYJ(lun,'NAXIS2',NCAfe,comment,ierr)
                         DO  i = 1 , NCAfe
                                CALL FTGCVJ(lun,1,i,1,1,0,IDNr(i),qanyf,ierr)
                                CALL FTGCVE(lun,2,i,1,1,0.,CDRcafes(i),qanyf,ierr)
                         ENDDO
                         CALL FTCLOS(lun,ierr)
                         IF ( ierr .NE. 0 ) THEN
                            if( myID == 0 ) WRITE(*,*) 'Failed to get data from mekal6.dat' 
                            STOP
                         ENDIF 
	                 DO i=1,L_CAFE
	                 cdrcafe(i)=cdrcafes(i)
	                 ENDDO
                    ENDIF
		  
!-----------------------------------------------------
            ! Second Read/Calculate telescope ARF and RMF function 


                    OPEN(unit=5,file=filARF,status='old')
                    OPEN(unit=15,file=filRMF,status='old')
                    OPEN(unit=20,file=filBG,status='old')      
                    OPEN(unit=25,file=filevent,status='old')  
		
                    DO j=1,xrayNbin
                         READ(5,*) ARF(j)
                   ENDDO  
		
           
                   DO j=1,xrayNbin
                        DO i=1,xrayNch 
                             READ(15,*) Rmf(j,i) 
			     TRM(j,i) =  Rmf(j,i)*Arf(j)	
                        ENDDO
                   ENDDO

                    DO i=1, LENx
                          READ(20,*)xrayBG_obs(i)      
                          READ(25,*)xrayCobs(i)
			
                    ENDDO     
	       
        
                 RETURN

       END SUBROUTINE ReadInData

END MODULE ReadWrite
