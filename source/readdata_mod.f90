module readdata_mod
use common_mod

!subroutines to read in data files to arrays at beginning of simulation

contains

    subroutine readData()

        implicit none
        integer :: i, ios                    ! counter, io state
        integer :: iup, ilow                 ! H counters
        integer :: elem_in, ion_in, ncore_in, nfit_in, k
        integer :: env_stat
        character(len=32) :: dr_env, rr_env
        real, dimension(9) :: c_in, e_in
        real :: a_rr_in, b_rr_in, t0_rr_in, t1_rr_in, c_rr_in, t2_rr_in

    ! initialise all the arrays that we are reading in

        y_dat = 0.
        Ay_dat = 0.
        HIRecLineData = 0.
        HeIILymanData = 0.
        HeIILymanNu = 0.
        HeIIRecLinedata = 0.
        direc_coeffs%elem = 0
        direc_coeffs%n = 0
        direc_coeffs%a = 0.
        direc_coeffs%b = 0.
        direc_coeffs%c = 0.
        direc_coeffs%d = 0.
        direc_coeffs%f = 0.
        direc_coeffs%g = 0
        aldropequi_coeffs%elem = 0
        aldropequi_coeffs%n = 0
        aldropequi_coeffs%a = 0
        aldropequi_coeffs%b = 0
        aldropequi_coeffs%t0 = 0.
        aldropequi_coeffs%t1 = 0.
        do elem_in = 1, nElements
           do ion_in = 1, nElements
              badnell_dr_coeffs(elem_in, ion_in)%nfit = 0
              badnell_dr_coeffs(elem_in, ion_in)%c = 0.
              badnell_dr_coeffs(elem_in, ion_in)%e = 0.
              badnell_rr_coeffs(elem_in, ion_in)%a = 0.
              badnell_rr_coeffs(elem_in, ion_in)%b = 0.
              badnell_rr_coeffs(elem_in, ion_in)%t0 = 0.
              badnell_rr_coeffs(elem_in, ion_in)%t1 = 0.
              badnell_rr_coeffs(elem_in, ion_in)%c = 0.
              badnell_rr_coeffs(elem_in, ion_in)%t2 = 0.
              badnell_rr_coeffs(elem_in, ion_in)%defined = .false.
           end do
        end do
        lgBadnellLoaded = .false.
        lgBadnellRRLoaded = .false.

    ! read in rates from data/HeI2phot.dat

        open(unit = 93,  action="read", file = PREFIX//"/share/mocassin/data/HeI2phot.dat", status = "old", position = "rewind", iostat=ios)
        if (ios /= 0) then
            print*, "! readData: can't open file: ",PREFIX,"/share/mocassin/data/HeI2phot.dat"
            stop
        end if
        do i = 1, 41
            read(unit = 93, fmt = *) y_dat(i), Ay_dat(i)
        end do

        close(93)

        ! read in HI recombination lines [e-25 ergs*cm^3/s]
        ! (Storey and Hummer MNRAS 272(1995)41)

        open(unit = 94,  action="read", file = PREFIX//"/share/mocassin/data/r1b0100.dat", status = "old", position = "rewind", iostat=ios)
        if (ios /= 0) then
            print*, "! readData: can't open file: ",PREFIX,"/share/mocassin/data/r1b0100.dat"
            stop
        end if

        do iup = 15, 3, -1
            read(94, fmt=*) (HIRecLinedata(iup, ilow), ilow = 2, min(8, iup-1))
        end do

        close(94)

        ! read in HeII Lyman line ratios up to level n=5 [e-25 ergs*cm^3/s]
        ! (Storey and Hummer MNRAS 272(1995)41)

        open(unit = 98,  action="read", file = PREFIX//"/share/mocassin/data/r2a0100.dat", status = "old", position = "rewind", iostat=ios)
        if (ios /= 0) then
            print*, "! readData: can't open file: ",PREFIX,"/share/mocassin/data/r2a0100.dat"
            stop
        end if

        do i = 1, NHeIILyman
            read(98, fmt=*) HeIILymanData(i), HeIILymanNu(i)
        end do

        close(98)

        ! read in HeII recombination lines [e-25 ergs*cm^3/s]
        ! (Storey and Hummer MNRAS 272(1995)41)

        open(unit = 95,  action="read", file = PREFIX//"/share/mocassin/data/r2b0100.dat", status = "old", position = "rewind", iostat=ios)
        if (ios /= 0) then
            print*, "! readData: can't open file: ",PREFIX,"/share/mocassin/data/r2b0100.dat"
            stop
        end if
        do iup = 30, 3, -1
            read(95, fmt=*) (HeIIRecLinedata(iup, ilow), ilow = 2, min(16, iup-1))
        end do

        close(95)

        ! dielectronic recombination coefficients
        call get_environment_variable("MOCASSIN_DR", dr_env, status=env_stat)
        if (env_stat == 0 .and. (trim(dr_env) == "legacy" .or. trim(dr_env) == "old")) then
           ios = -1
        else
           open (unit=18, file=PREFIX//'/share/mocassin/data/badnell_dr.dat', status='old', position='rewind', iostat = ios, action="read")
        end if

        if (ios == 0) then
           do
              read(unit=18, fmt=*, iostat=ios) elem_in, ion_in, ncore_in, nfit_in, (c_in(k), k=1,9), (e_in(k), k=1,9)
              if (ios /= 0) exit
              if (elem_in >= 1 .and. elem_in <= nElements .and. ion_in >= 1 .and. ion_in <= nElements) then
                 badnell_dr_coeffs(elem_in, ion_in)%nfit = nfit_in
                 badnell_dr_coeffs(elem_in, ion_in)%c = c_in
                 badnell_dr_coeffs(elem_in, ion_in)%e = e_in
                 lgBadnellLoaded = .true.
              end if
           end do
           close(18)
           print*, "! readData: using Badnell dielectronic recombination dataset (Cloudy c25.00)"
        else
           ! Fall back to legacy dielectronic.dat (Nussbaumer & Storey)
           open (unit=18, file=PREFIX//'/share/mocassin/data/dielectronic.dat', status='old',position='rewind', iostat = ios, action="read")
           do i = 1, 25
              read(unit=18, fmt=*, iostat=ios) direc_coeffs(i)%elem, direc_coeffs(i)%n, direc_coeffs(i)%a, direc_coeffs(i)%b, direc_coeffs(i)%c, direc_coeffs(i)%d, direc_coeffs(i)%f, direc_coeffs(i)%g
              if (ios < 0) exit ! end of file reached
           enddo
           close(18)
           print*, "! readData: using legacy dielectronic recombination dataset (Nussbaumer & Storey 1983-86)"
        end if

        ! high temperature dielectronic recombination coefficients from
        ! Aldrovandi and Pequignot 1973

        open (unit=17, file=PREFIX//'/share/mocassin/data/aldrovandi.dat', status='old',position='rewind', iostat = ios, action="read")
        if (ios /= 0) then
           print*, "! readData: can't open file ",PREFIX,"/share/mocassin/data/alrovandi.dat"
           stop
        end if

        do i = 1, 167
            read(unit=17, fmt=*, iostat=ios) aldropequi_coeffs(i)%elem, aldropequi_coeffs(i)%n, aldropequi_coeffs(i)%a, aldropequi_coeffs(i)%b, aldropequi_coeffs(i)%t0, aldropequi_coeffs(i)%t1
            if (ios < 0) exit ! end of file reached
        end do

        close(17)

        ! radiative recombination coefficients from Badnell (Cloudy c25.00)
        call get_environment_variable("MOCASSIN_RR", rr_env, status=env_stat)
        if (env_stat == 0 .and. (trim(rr_env) == "legacy" .or. trim(rr_env) == "old")) then
           ios = -1
        else
           open (unit=19, file=PREFIX//'/share/mocassin/data/badnell_rr.dat', status='old', position='rewind', iostat = ios, action="read")
        end if

        if (ios == 0) then
           do
              read(unit=19, fmt=*, iostat=ios) elem_in, ion_in, ncore_in, a_rr_in, b_rr_in, t0_rr_in, t1_rr_in, c_rr_in, t2_rr_in
              if (ios /= 0) exit
              if (elem_in >= 1 .and. elem_in <= nElements .and. ion_in >= 1 .and. ion_in <= nElements) then
                 badnell_rr_coeffs(elem_in, ion_in)%a = a_rr_in
                 badnell_rr_coeffs(elem_in, ion_in)%b = b_rr_in
                 badnell_rr_coeffs(elem_in, ion_in)%t0 = t0_rr_in
                 badnell_rr_coeffs(elem_in, ion_in)%t1 = t1_rr_in
                 badnell_rr_coeffs(elem_in, ion_in)%c = c_rr_in
                 badnell_rr_coeffs(elem_in, ion_in)%t2 = t2_rr_in
                 badnell_rr_coeffs(elem_in, ion_in)%defined = .true.
                 lgBadnellRRLoaded = .true.
              end if
           end do
           close(19)
           print*, "! readData: using Badnell radiative recombination dataset (Cloudy c25.00)"
        else
           print*, "! readData: using legacy radiative recombination dataset (Verner & Ferland 1996)"
        end if

    end subroutine readData

end module readdata_mod

