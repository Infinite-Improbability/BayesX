module ReadInputs

   use params

   implicit none

contains

   integer function getinputs(infilename)

      implicit none

      integer ios, funit
      character*500 infilename, tag

      getinputs = 0        !error flag

      funit = 22
      open (unit=funit, file=infilename, status='old')
      do
         read (funit, *, IOSTAT=ios) tag
         if (ios < 0) exit

         tag = trim(tag)
         if (scan(tag, '#', .false.) == 0) then
            cycle
         elseif (tag == '#filion') then
            read (funit, *) filion
            filion = trim(filion)
         elseif (tag == '#filrec') then
            read (funit, *) filrec
            filrec = trim(filrec)
         elseif (tag == '#filkar') then
            read (funit, *) filkar
            filkar = trim(filkar)
         elseif (tag == '#filcon') then
            read (funit, *) filcon
            filcon = trim(filcon)
         elseif (tag == '#fillin') then
            read (funit, *) fillin
            fillin = trim(fillin)
         elseif (tag == '#filcaf') then
            read (funit, *) filcaf
            filcaf = trim(filcaf)
         elseif (tag == '#filARF') then
            read (funit, *) filARF
            filARF = trim(filARF)
         elseif (tag == '#filRMF') then
            read (funit, *) filRMF
            filRMF = trim(filRMF)
         elseif (tag == '#filBG') then
            read (funit, *) filBG
            filBG = trim(filBG)
         elseif (tag == '#filevent') then
            read (funit, *) filevent
            filevent = trim(filevent)
         elseif (tag == '#filmask') then
            read (funit, *) filmask
            filmask = trim(filmask)
         elseif (tag == '#XrayTelescope') then
            read (funit, *) XrayTelescope
            XrayTelescope = trim(XrayTelescope)
            if (len(XrayTelescope) > 3) then
               if (myID == 0) write (*, *) 'ERROR: Telescope name can not exceed 3 characters.'
               getinputs = 1
               exit
            end if
         elseif (tag == '#n') then
            read (funit, *) n
         elseif (tag == '#nx') then
            read (funit, *) xraynx
            xraynxcentre = xraynx/2 + 1
         elseif (tag == '#ny') then
            read (funit, *) xrayny
            xraynycentre = xrayny/2 + 1
         elseif (tag == '#Aeffave') then
            read (funit, *) Aeffave
         elseif (tag == '#xraycell') then
            read (funit, *) xraycell
         elseif (tag == '#xrayNbin') then
            read (funit, *) xrayNbin
         elseif (tag == '#xrayNch') then
            read (funit, *) xrayNch
         elseif (tag == '#xrayEmin') then
            read (funit, *) xrayEmin
         elseif (tag == '#xrayEmax') then
            read (funit, *) xrayEmax
         elseif (tag == '#sexpotime') then
            read (funit, *) sexpotime
         elseif (tag == '#bexpotime') then
            read (funit, *) bexpotime
         elseif (tag == '#NHcol') then
            read (funit, *) N_H_col
         elseif (tag == '#xrayBG_model') then
            read (funit, *) xrayBG_predmax
         elseif (tag == '#IS') then
            read (funit, *) n_IS
         elseif (tag == '#multimodal') then
            read (funit, *) n_mmodal
         elseif (tag == '#nlive') then
            read (funit, *) n_nlive
         elseif (tag == '#eff') then
            read (funit, *) n_efr
         elseif (tag == '#tol') then
            read (funit, *) n_tol
         elseif (tag == '#updint') then
            read (funit, *) n_updint
         elseif (tag == '#maxmodes') then
            read (funit, *) n_maxModes
         elseif (tag == '#nCdims') then
            read (funit, *) clusterDims
         elseif (tag == '#seed') then
            read (funit, *) n_rseed
         elseif (tag == '#root') then
            read (funit, *) n_root
            n_root = trim(n_root)
         elseif (tag == '#cluster_model') then
            read (funit, *) GasModel
         elseif (tag == '#x_prior') then
            read (funit, *) Geo_PriorType(1, 1), Geo_Prior(1, 1, 1), Geo_Prior(1, 1, 2)
         elseif (tag == '#y_prior') then
            read (funit, *) Geo_PriorType(1, 2), Geo_Prior(1, 2, 1), Geo_Prior(1, 2, 2)
         elseif (tag == '#m200_prior') then
            read (funit, *) Gas_PriorType(1, 1), Gas_Prior(1, 1, 1), Gas_Prior(1, 1, 2)
         elseif ((tag == '#fgas200_prior') .and. (GasModel .ne. 3)) then
            read (funit, *) Gas_PriorType(1, 2), Gas_Prior(1, 2, 1), Gas_Prior(1, 2, 2)
         elseif ((tag == '#a_GNFW_prior') .and. (GasModel .ne. 3)) then
            read (funit, *) Gas_PriorType(1, 3), Gas_Prior(1, 3, 1), Gas_Prior(1, 3, 2)
         elseif ((tag == '#b_GNFW_prior') .and. (GasModel .ne. 3)) then
            read (funit, *) Gas_PriorType(1, 4), Gas_Prior(1, 4, 1), Gas_Prior(1, 4, 2)
         elseif ((tag == '#c_GNFW_prior') .and. (GasModel .ne. 3)) then
            read (funit, *) Gas_PriorType(1, 5), Gas_Prior(1, 5, 1), Gas_Prior(1, 5, 2)
         elseif ((tag == '#c500_GNFW_prior') .and. (GasModel .ne. 3)) then
            read (funit, *) Gas_PriorType(1, 6), Gas_Prior(1, 6, 1), Gas_Prior(1, 6, 2)
         elseif ((tag == '#alpha_model2_prior') .and. (GasModel == 2)) then
            read (funit, *) Gas_PriorType(1, 7), Gas_Prior(1, 7, 1), Gas_Prior(1, 7, 2)
         elseif ((tag == '#gamma0_poly_prior') .and. (GasModel == 3)) then
            read (funit, *) Gas_PriorType(1, 8), Gas_Prior(1, 8, 1), Gas_Prior(1, 8, 2)
         elseif ((tag == '#gammaR_poly_prior') .and. (GasModel == 3)) then
            read (funit, *) Gas_PriorType(1, 9), Gas_Prior(1, 9, 1), Gas_Prior(1, 9, 2)
         elseif ((tag == '#t0_poly_prior') .and. (GasModel == 3)) then
            read (funit, *) Gas_PriorType(1, 10), Gas_Prior(1, 10, 1), Gas_Prior(1, 10, 2)
         elseif (tag == '#z_Prior') then
            read (funit, *) z_PriorType(1), z_Prior(1, 1), z_Prior(1, 2)
         elseif (tag == '#mass_function') then
            read (funit, *) mass_function
         elseif (tag == '#rauto') then
            read (funit, *) rauto
         elseif (tag == '#rmin') then
            read (funit, *) rmin
         elseif (tag == '#rmax') then
            read (funit, *) rmax
         elseif (tag == '#rlimit') then
            read (funit, *) rlimit
         end if
      end do

      close (funit)

      !check if all the compulsory inputs have been defined
      if (filion == '') then
         getinputs = 1
         if (myID == 0) write (*, *) "ERROR: filion not set in the input file"
      end if
      if (filrec == '') then
         getinputs = 1
         if (myID == 0) write (*, *) "ERROR: filrec not set in the input file"
      end if
      if (filkar == '') then
         getinputs = 1
         if (myID == 0) write (*, *) "ERROR: filkar not set in the input file"
      end if
      if (filcon == '') then
         getinputs = 1
         if (myID == 0) write (*, *) "ERROR: filcon not set in the input file"
      end if
      if (fillin == '') then
         getinputs = 1
         if (myID == 0) write (*, *) "ERROR: fillin not set in the input file"
      end if
      if (filcaf == '') then
         getinputs = 1
         if (myID == 0) write (*, *) "ERROR: filcaf not set in the input file"
      end if
      if (filARF == '') then
         getinputs = 1
         if (myID == 0) write (*, *) "ERROR: filARF not set in the input file"
      end if
      if (filRMF == '') then
         getinputs = 1
         if (myID == 0) write (*, *) "ERROR: filRMF not set in the input file"
      end if
      if (filBG == '') then
         getinputs = 1
         if (myID == 0) write (*, *) "ERROR: filBG not set in the input file"
      end if
      if (filevent == '') then
         getinputs = 1
         if (myID == 0) write (*, *) "ERROR: filevent not set in the input file"
      end if
      if (n_root == '') then
         getinputs = 1
         if (myID == 0) write (*, *) "ERROR: root not set in the input file"
      end if
      if (xraynx == 0) then
         getinputs = 1
         if (myID == 0) write (*, *) "ERROR: nx not set in the input file"
      end if
      if (xrayny == 0) then
         getinputs = 1
         if (myID == 0) write (*, *) "ERROR: ny not set in the input file"
      end if
      if (xrayNbin == 0) then
         getinputs = 1
         if (myID == 0) write (*, *) "ERROR: xrayNbin not set in the input file"
      end if
      if (xrayNch == 0) then
         getinputs = 1
         if (myID == 0) write (*, *) "ERROR: xrayNch not set in the input file"
      end if

   end function getinputs

end module ReadInputs
