!=====================================================================
!
!                            BAYES-X
!
!  Bayesian inference tool for the analysis of X-ray observations of galaxy clusters
!
!  Author:        M. Olamaie et al
!                 September 2014
!
!=======================================================================

PROGRAM BayesX

   USE params
   USE like
   USE ReadWrite
   USE massfunction
   USE nestwrapper
   USE Nested
   USE ReadInputs

   IMPLICIT NONE

   INTEGER                          :: index, i, j, k, context, eflag
   REAL*8                           :: Evidence, ev, temp(4), pos(2), slike, s_nullev, d1, d2
   LOGICAL                                :: flag
   CHARACTER(len=500)        :: infile

!-----------------------------------------------------------------------

#ifdef MPI
   !MPI initializations
   CALL MPI_INIT(errcode)
   IF (errcode /= MPI_SUCCESS) THEN
      WRITE (*, *) 'Error starting MPI program. Terminating.'
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode)
   END IF
   call MPI_COMM_RANK(MPI_COMM_WORLD, myID, errcode)
#endif

!-----------------------------------------------------------------------

   eflag = 0
   if (myID == 0) then
      if (iargc() /= 1) then
         write (*, *) "Usage: ./BayesX [input file]"
         eflag = 1
      else
         call getarg(1, infile)
         inquire (file=infile, exist=flag)
         if (.not. flag) then
            write (*, *) "ERROR: File ", trim(infile), " does not exist."
            eflag = 1
         end if
      end if
   end if

#ifdef MPI
   call MPI_BARRIER(MPI_COMM_WORLD, errcode)
   call MPI_BCAST(eflag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, errcode)
   call MPI_BCAST(infile, 500, MPI_CHARACTER, 0, MPI_COMM_WORLD, errcode)
#endif

   if (eflag == 0) then
      eflag = getinputs(infile)
   end if
   !write (*, *) filion

   if (eflag == 0) then
      LENx = xraynx*xrayny*xrayNch
      xrayDeltaE = (xrayEmax - xrayEmin)/xrayNbin

      allocate (photar(xrayNbin), xrayE1(xrayNbin), xrayE2(xrayNbin), xrayFluxCoeff(xrayNbin), xrayFlux(xrayNbin), &
                xraypredFlux(xrayNch), Arf(xrayNbin), Rmf(xrayNbin, xrayNch), TRM(xrayNbin, xrayNch), xrayBG(LENx), &
                xrayCobs(LENx), xrayBG_obs(LENX), xrayCmap(xrayNch, xraynx, xrayny), xrayCpred(LENx), r(n), T(n), Rhogas(n), logr(n))

      index = 1
      if (myID == 0) CALL Welcome
      CALL ReadInData
      WRITE (*, *)
      eflag = Initialise(index)
      if (eflag == 0) then
         IF (edim > 15) n_pWrap = 1
         CALL srand(n_rseed)

         call nestRun(n_IS, n_mmodal, n_ceff, n_nlive, n_tol, n_efr, edim, n_totPar, 2, n_maxModes, n_updInt, -1d99, &
                      n_root, n_rseed, n_pWrap, n_fb, .true., .true., .false., -1d10, -1, getloglike, dumper, context)
      end if
   end if

!-----------------------------------------------------------------------

#ifdef MPI
   CALL MPI_FINALIZE(errcode)
#endif

END PROGRAM BayesX
