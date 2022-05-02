MODULE ReadWrite
   USE params
   USE utilities
   USE constants
   USE MassModels

CONTAINS

   SUBROUTINE ReadInData

      IMPLICIT NONE

      INTEGER                         ::  i, j, m, iend
      INTEGER                         ::  ii, jj, mm
      REAL*8                          ::  ra, dec, ra2, dec2, rao, deco
      CHARACTER(LEN=100)    :: string

      ! Reading X ray data

      ! First read data for MEKAL Model

      IF (init .EQ. 0) THEN
         init = 1
         CALL ftopen(lun, filion, 0, indexy, ierr)
         CALL FTMRHD(lun, 1, indexy, ierr)
         DO i = 1, NREC
            CALL FTGCVE(lun, 1, i, 1, 256, 0., AIOns((i - 1)*256 + 1), qanyf, ierr)
         END DO
         CALL FTCLOS(lun, ierr)
         IF (ierr .NE. 0) THEN
            if (myID == 0) WRITE (*, *) 'ERROR: Failed to get data from mekal1.dat'
            STOP
         END IF
         DO i = 1, NAMAX
            aion(i) = AIOns(i)
         END DO
         CALL ftopen(lun, filrec, 0, indexy, ierr)
         CALL FTMRHD(lun, 1, indexy, ierr)
         DO i = 1, NUMION
            CALL FTGCVE(lun, 1, i, 1, 13, 0., AREcs(1, i), qanyf, ierr)
         END DO
         CALL FTCLOS(lun, ierr)
         IF (ierr .NE. 0) THEN
            if (myID == 0) WRITE (*, *) 'Failed to get data from mekal2.dat'
            STOP
         END IF
         DO i = 1, 13
            DO j = 1, NUMION
               arec(i, j) = AREcs(i, j)
            END DO
         END DO
         CALL ftopen(lun, filkar, 0, indexy, ierr)
         CALL FTMRHD(lun, 1, indexy, ierr)
         DO i = 1, NUU
            CALL FTGCVE(lun, 1, i, 1, NG, 0., GAsing(1, i), qanyf, ierr)
         END DO
         CALL FTCLOS(lun, ierr)
         IF (ierr .NE. 0) THEN
            if (myID == 0) WRITE (*, *) 'Failed to get data from mekal3.dat'
            STOP
         END IF
         DO i = 1, NG
            DO j = 1, NUU
               ga(i, j) = gAsing(i, j)
            END DO
         END DO
         CALL ftopen(lun, filcon, 0, indexy, ierr)
         CALL FTMRHD(lun, 1, indexy, ierr)
         DO i = 1, NUMION
            CALL FTGCVE(lun, 1, i, 1, 6, 0., Psing(1, i), qanyf, ierr)
         END DO
         CALL FTCLOS(lun, ierr)
         IF (ierr .NE. 0) THEN
            if (myID == 0) WRITE (*, *) 'Failed to get data from mekal4.dat'
            STOP
         END IF
         DO i = 1, 6
            DO j = 1, NUMION
               p(i, j) = Psing(i, j)
            END DO
         END DO
         CALL ftopen(lun, fillin, 0, indexy, ierr)
         CALL FTMRHD(lun, 1, indexy, ierr)
         CALL FTGKYJ(lun, 'NAXIS2', NE, comment, ierr)
         DO i = 1, NE
            CALL FTGCVI(lun, 1, i, 1, 1, 0, NRL(i), qanyf, ierr)
            CALL FTGCVI(lun, 2, i, 1, 9, 0, LP(1, i), qanyf, ierr)
            CALL FTGCVS(lun, 3, i, 1, 1, ' ', TRAns(i), qanyf, ierr)
            CALL FTGCVE(lun, 4, i, 1, 14, 0., RPsing(1, i), qanyf, ierr)
         END DO
         CALL FTCLOS(lun, ierr)
         IF (ierr .NE. 0) THEN
            if (myID == 0) WRITE (*, *) 'Failed to get data from mekal5.dat'
            STOP
         END IF
         DO i = 1, 14
            DO j = 1, NSTO
               rp(i, j) = RPsing(i, j)
            END DO
         END DO
         CALL ftopen(lun, filcaf, 0, indexy, ierr)
         CALL FTMRHD(lun, 1, indexy, ierr)
         CALL FTGKYJ(lun, 'NAXIS2', NCAfe, comment, ierr)
         DO i = 1, NCAfe
            CALL FTGCVJ(lun, 1, i, 1, 1, 0, IDNr(i), qanyf, ierr)
            CALL FTGCVE(lun, 2, i, 1, 1, 0., CDRcafes(i), qanyf, ierr)
         END DO
         CALL FTCLOS(lun, ierr)
         IF (ierr .NE. 0) THEN
            if (myID == 0) WRITE (*, *) 'Failed to get data from mekal6.dat'
            STOP
         END IF
         DO i = 1, L_CAFE
            cdrcafe(i) = cdrcafes(i)
         END DO
      END IF

!-----------------------------------------------------
      ! Second Read/Calculate telescope ARF and RMF function

      OPEN (unit=5, file=filARF, status='old')
      OPEN (unit=15, file=filRMF, status='old')
      OPEN (unit=20, file=filBG, status='old')
      OPEN (unit=25, file=filevent, status='old')

      DO j = 1, xrayNbin
         READ (5, *) ARF(j)
      END DO

      DO j = 1, xrayNbin
         DO i = 1, xrayNch
            READ (15, *) Rmf(j, i)
            TRM(j, i) = Rmf(j, i)*Arf(j)
         END DO
      END DO

      DO i = 1, LENx
         READ (20, *) xrayBG_obs(i)
         READ (25, *) xrayCobs(i)

      END DO

      RETURN

   END SUBROUTINE ReadInData

!=======================================================================
   subroutine write_paramnames
      implicit none

      integer ival, i, j, ii, ival1
      character(len=3) :: js, iis
      character(len=4) :: ext
      character(len=50) :: parname,fmt,parunit
      integer, parameter :: params_unit = 21

      open(unit=params_unit, form='formatted', file=trim(n_root)//'.paramnames', status='replace')
      ival1 = 0
      do ival = 1, Ndim

            j=(ival-1)/NPars+1
            write(js,'(I0)') j
            i=ival-(j-1)*NPars
            parunit=''
            if(i <= NGeoPars) then
               ii = i
               if (Geo_PriorType(j,ii)<=0) cycle
               if (ii == 1) then
                  parname='x_0'
                  if(NAtoms.gt.1) parname='x_{0,'//trim(js)//'}'
                  parunit='\mathrm{arcsec}'
               elseif (ii == 2) then
                  parname='y_0'
                  if(NAtoms.gt.1) parname='y_{0,'//trim(js)//'}'
                  parunit='\mathrm{arcsec}'
               elseif (ii == 3) then
                  parname='\theta'
                  if(NAtoms.gt.1) parname=trim(parname)//'_{'//trim(js)//'}'
                  parunit='\mathrm{deg}'
               elseif (ii == 4) then
                  parname='f'
                  if(NAtoms.gt.1) parname=trim(parname)//'_{'//trim(js)//'}'
               endif
            elseif (i <= NGeoPars+1 ) then
               if (z_PriorType(j)<=0) cycle
               parname='z'
               if(NAtoms.gt.1) parname=trim(parname)//'_{'//trim(js)//'}'
            else if( Gas == 1 .and. i <= ( NGeoPars + 1 + NGasPars * Gas ) ) then
               ii = i - NGeoPars - 1
               if (Gas_PriorType(j,ii)<=0) cycle
               if (ii==1) then
                   parname='M_{T,200}'
                   if(NAtoms.gt.1) parname='M_{T,200,'//trim(js)//'}'
                   parunit='M_{\odot}'
               elseif (ii==2) then
                   parname='f_{\mathrm{gas},200}'
                   if(NAtoms.gt.1) parname='f_{\mathrm{gas},200,'//trim(js)//'}'
               elseif (ii==3) then
                   parname='\gamma'
                   if(NAtoms.gt.1) parname=trim(parname)//'_{'//trim(js)//'}'
               elseif (ii==4) then
                   parname='\alpha'
                   if(NAtoms.gt.1) parname=trim(parname)//'_{'//trim(js)//'}'
               elseif (ii==5) then
                   parname='\beta'
                   if(NAtoms.gt.1) parname=trim(parname)//'_{'//trim(js)//'}'
               elseif (ii==6) then
                   parname='c_{500}'
                   if(NAtoms.gt.1) parname='c_{500,'//trim(js)//'}'
	       elseif (ii==7) then
                   parname='\alpha_{Ein}'
                   if(NAtoms.gt.1) parname='\alpha_{Ein,'//trim(js)//'}'
               endif
            endif
            ival1 = ival1 + 1
            fmt='(2X,A1,I3.3,6X,A)'
            if (len(trim(parunit)).gt.0) parname=trim(parname)//' / '//trim(parunit)
            write(params_unit,fmt) 'p', ival1, trim(parname)
      enddo

      ! Auxiliary parameters - name ending in '*' means getdist.py interprets them as derived
      ext='}'
      do j = 1, NAtoms
        write(js,'(I0)') j
        if (NAtoms > 1) then
          ext=','//trim(js)//'}'
        endif
        fmt='(2X,A1,I3.3,A1,5X,A)'
        if(NAtoms==1) then
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+1, '*', 'D / \mathrm{Mpc}'
        else
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+1, '*', 'D_{'//trim(js)//'} / \mathrm{Mpc}'
        endif
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+2, '*', 'r_{s'//trim(ext)//' / \mathrm{Mpc}'
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+3, '*', '\rho_{s'//trim(ext)
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+4, '*', 'r_{p'//trim(ext)//' / \mathrm{Mpc}'
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+5, '*', 'P_{ei'//trim(ext)
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+6, '*', 'r_{2500'//trim(ext)//' / \mathrm{Mpc}'
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+7, '*', 'c_{2500'//trim(ext)
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+8, '*', 'M_{g,2500'//trim(ext)//' / M_{\odot}'
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+9, '*', 'M_{T,2500'//trim(ext)//' / M_{\odot}'
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+10, '*', 'f_{g,2500'//trim(ext)
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+11, '*', 'T_{g,2500'//trim(ext)	 	 	 		
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+12, '*', 'n_{e,2500'//trim(ext)	 	 	 		
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+13, '*', 'K_{e,2500'//trim(ext)	 	 	 		
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+14, '*', 'P_{e,2500'//trim(ext)	 	 	 		
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+15, '*', 'r_{500'//trim(ext)//' / \mathrm{Mpc}'
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+16, '*', 'c_{500'//trim(ext)
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+17, '*', 'M_{g,500'//trim(ext)//' / M_{\odot}'
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+18, '*', 'M_{T,500'//trim(ext)//' / M_{\odot}'
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+19, '*', 'f_{g,500'//trim(ext)
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+20, '*', 'T_{g,500'//trim(ext)	 	 	 		
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+21, '*', 'n_{e,500'//trim(ext)	 	 	 		
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+22, '*', 'K_{e,500'//trim(ext)	 	 	 		
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+23, '*', 'P_{e,500'//trim(ext)	 	 	 		
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+24, '*', 'r_{200'//trim(ext)//' / \mathrm{Mpc}'
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+25, '*', 'c_{200'//trim(ext)
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+26, '*', 'M_{g,200'//trim(ext)//' / M_{\odot}'
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+27, '*', 'M_{T,200a'//trim(ext)//' / M_{\odot}' ! Just give it a different label from the sampling parameter to avoid confusing getdist
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+28, '*', 'f_{g,200'//trim(ext)
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+29, '*', 'T_{g,200'//trim(ext)	 	 	 		
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+30, '*', 'n_{e,200'//trim(ext)	 	 	 		
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+31, '*', 'K_{e,200'//trim(ext)	 	 	 		
        write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+32, '*', 'P_{e,200'//trim(ext)	 	 	 		
      enddo  

      close(params_unit)
    end subroutine write_paramnames
!=======================================================================

END MODULE ReadWrite
