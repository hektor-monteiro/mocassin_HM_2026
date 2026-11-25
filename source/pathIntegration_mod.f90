! Copyright (C) 2005 Barbara Ercolano
!
! Version 2.02
module pathIntegration_mod
    use common_mod
    use interpolation_mod
    use vector_mod

    private :: initPathTauStart, stepPathTauCell

    contains

    ! Locate the starting cell indices (xP, yP, zP) for a ray originating at aVec
    ! and validate that the starting position lies inside the single-grid domain.
    subroutine initPathTauStart(grid, aVec, xP, yP, zP, callerName)
        implicit none

        type(grid_type),  intent(in)  :: grid(*)
        type(vector),     intent(in)  :: aVec
        integer,          intent(out) :: xP, yP, zP
        character(len=*), intent(in)  :: callerName

        real :: R_cyl

        if (lg2D .and. lgCylinder2D) then
           R_cyl = real(sqrt(dble(aVec%x)**2 + dble(aVec%y)**2))
           call locate(grid(1)%xAxis, R_cyl, xP)
           if (xP == 0) xP = 1
           if (xP < grid(1)%nx) then
              if (R_cyl >= (grid(1)%xAxis(xP) + grid(1)%xAxis(xP+1))/2.) xP = xP + 1
           end if

           yP = 1

           call locate(grid(1)%zAxis, aVec%z, zP)
           if (zP == 0 .and. lgSymmetricXYZ .and. aVec%z >= grid(1)%zAxis(1)) zP = 1
           if (zP < grid(1)%nz) then
              if (aVec%z >= (grid(1)%zAxis(zP) + grid(1)%zAxis(zP+1))/2.) zP = zP + 1
           end if

           if ((xP <= 0) .or. (xP > grid(1)%nx)) then
              print*, trim(callerName), ": starting position in R is outside the grid", xP, yP, zP
              stop
           end if
           if ((zP <= 0) .or. (zP > grid(1)%nz)) then
              print*, trim(callerName), ": starting position in z is outside the grid", xP, yP, zP
              stop
           end if
           return
        end if

        ! Locate the starting position in 3D or planar 2D Cartesian grid
        if (lg2D .and. perpAxis2D == "X") then
           xP = 1
        else
           call locate(grid(1)%xAxis, aVec%x, xP)
           if (xP < grid(1)%nx) then
              if (aVec%x >= (grid(1)%xAxis(xP) + grid(1)%xAxis(xP+1))/2.) xP = xP + 1
           end if
        end if

        if (lg2D .and. perpAxis2D == "Y") then
           yP = 1
        else
           call locate(grid(1)%yAxis, aVec%y, yP)
           if (yP < grid(1)%ny) then
              if (aVec%y >= (grid(1)%yAxis(yP) + grid(1)%yAxis(yP+1))/2.) yP = yP + 1
           end if
        end if

        if (lg2D .and. perpAxis2D == "Z") then
           zP = 1
        else
           call locate(grid(1)%zAxis, aVec%z, zP)
           if (zP < grid(1)%nz) then
              if (aVec%z >= (grid(1)%zAxis(zP) + grid(1)%zAxis(zP+1))/2.) zP = zP + 1
           end if
        end if

        if (.not. (lg2D .and. perpAxis2D == "X")) then
           if ((xP <= 0) .or. (xP > grid(1)%nx)) then
              print*, trim(callerName), ": starting position in x is outside the grid", xP, yP, zP
              stop
           end if
        end if
        if (.not. (lg2D .and. perpAxis2D == "Y")) then
           if ((yP <= 0) .or. (yP > grid(1)%ny)) then
              print*, trim(callerName), ": starting position in y is outside the grid", xP, yP, zP
              stop
           end if
        end if
        if (.not. (lg2D .and. perpAxis2D == "Z")) then
           if ((zP <= 0) .or. (zP > grid(1)%nz)) then
              print*, trim(callerName), ": starting position in z is outside the grid", xP, yP, zP
              stop
           end if
        end if

    end subroutine initPathTauStart


    ! Check outer/inner domain exits, apply physical symmetry reflections, and
    ! update cell indices (xP, yP, zP) for the current ray position rVec.
    subroutine stepPathTauCell(grid, rVec, vHat, xP, yP, zP, lgExit)
        implicit none

        type(grid_type), intent(in)    :: grid(*)
        type(vector),    intent(inout) :: rVec, vHat
        integer,         intent(inout) :: xP, yP, zP
        logical,         intent(out)   :: lgExit

        real    :: R_cyl
        logical :: lgHasInPlaneVel

        lgExit = .false.

        if (lg2D .and. lgCylinder2D) then
           R_cyl = real(sqrt(dble(rVec%x)**2 + dble(rVec%y)**2))

           if (R_cyl > grid(1)%xAxis(grid(1)%nx)) then
              lgExit = .true.
              return
           end if
           if (rVec%z > grid(1)%zAxis(grid(1)%nz)) then
              lgExit = .true.
              return
           end if
           if (.not. lgSymmetricXYZ .and. rVec%z < grid(1)%zAxis(1)) then
              lgExit = .true.
              return
           end if
           if (real(sqrt(dble(rVec%x)**2 + dble(rVec%y)**2 + dble(rVec%z)**2)) >= R_out .and. R_out > 0.) then
              lgExit = .true.
              return
           end if

           ! Equatorial mirror symmetry in true axisymmetric cylindrical (R, Z)
           if (lgSymmetricXYZ .and. rVec%z < grid(1)%zAxis(1)) then
              vHat%z = -vHat%z
              rVec%z = 2.0*grid(1)%zAxis(1) - rVec%z
           end if

           yP = 1
           call locate(grid(1)%xAxis, R_cyl, xP)
           if (xP == 0) xP = 1
           if (xP < grid(1)%nx) then
              if (R_cyl > (grid(1)%xAxis(xP) + grid(1)%xAxis(xP+1))/2.) xP = xP + 1
           else if (xP == grid(1)%nx) then
              if (R_cyl > grid(1)%xAxis(xP)) then
                 lgExit = .true.
                 return
              end if
           end if

           call locate(grid(1)%zAxis, rVec%z, zP)
           if (zP == 0) then
              if (lgSymmetricXYZ) then
                 zP = 1
              else
                 lgExit = .true.
                 return
              end if
           end if
           if (zP < grid(1)%nz) then
              if (rVec%z > (grid(1)%zAxis(zP) + grid(1)%zAxis(zP+1))/2.) zP = zP + 1
           else if (zP == grid(1)%nz) then
              if (rVec%z > grid(1)%zAxis(zP)) then
                 lgExit = .true.
                 return
              end if
           end if

           return
        end if

        ! Check if the path is still within the 3D or planar 2D domain
        if (.not. (lg2D .and. perpAxis2D == "X")) then
           if (rVec%x > grid(1)%xAxis(grid(1)%nx)) then
              lgExit = .true.
              return
           end if
        end if

        if (.not. (lg2D .and. perpAxis2D == "Y")) then
           if (rVec%y > grid(1)%yAxis(grid(1)%ny)) then
              lgExit = .true.
              return
           end if
        end if

        if (.not. (lg2D .and. perpAxis2D == "Z")) then
           if (rVec%z > grid(1)%zAxis(grid(1)%nz)) then
              lgExit = .true.
              return
           end if
        end if

        if (real(sqrt(dble(rVec%x)**2 + dble(rVec%y)**2 + dble(rVec%z)**2)) >= R_out .and. R_out > 0.) then
           lgExit = .true.
           return
        end if

        if (lgSymmetricXYZ) then
           if (.not. (lg2D .and. perpAxis2D == "X")) then
              if (rVec%x < grid(1)%xAxis(1)) then
                 vHat%x = -vHat%x
                 rVec%x = -rVec%x
                 call locate(grid(1)%xAxis, rVec%x, xP)
              end if
           end if
           if (.not. (lg2D .and. perpAxis2D == "Y")) then
              if (rVec%y < grid(1)%yAxis(1)) then
                 vHat%y = -vHat%y
                 rVec%y = -rVec%y
                 call locate(grid(1)%yAxis, rVec%y, yP)
              end if
           end if
           if (.not. (lg2D .and. perpAxis2D == "Z")) then
              if (rVec%z < grid(1)%zAxis(1)) then
                 vHat%z = -vHat%z
                 rVec%z = -rVec%z
                 call locate(grid(1)%zAxis, rVec%z, zP)
              end if
           end if
        end if

        ! Planar 2D perpendicular wall handling: never reflect rays with zero in-plane velocity
        if (lg2D) then
           select case (perpAxis2D)
           case ("Z")
              zP = 1
              lgHasInPlaneVel = (abs(vHat%x) > 1.e-20 .or. abs(vHat%y) > 1.e-20)
              if (rVec%z < minPerp2D) then
                 if ((lgSymmetricXYZ .or. lgSpecular2D) .and. lgHasInPlaneVel) then
                    vHat%z = -vHat%z
                    rVec%z = 2.0*minPerp2D - rVec%z
                 else
                    lgExit = .true.
                    return
                 end if
              else if (rVec%z > maxPerp2D) then
                 if (lgSpecular2D .and. lgHasInPlaneVel) then
                    vHat%z = -vHat%z
                    rVec%z = 2.0*maxPerp2D - rVec%z
                 else
                    lgExit = .true.
                    return
                 end if
              end if
           case ("Y")
              yP = 1
              lgHasInPlaneVel = (abs(vHat%x) > 1.e-20 .or. abs(vHat%z) > 1.e-20)
              if (rVec%y < minPerp2D) then
                 if ((lgSymmetricXYZ .or. lgSpecular2D) .and. lgHasInPlaneVel) then
                    vHat%y = -vHat%y
                    rVec%y = 2.0*minPerp2D - rVec%y
                 else
                    lgExit = .true.
                    return
                 end if
              else if (rVec%y > maxPerp2D) then
                 if (lgSpecular2D .and. lgHasInPlaneVel) then
                    vHat%y = -vHat%y
                    rVec%y = 2.0*maxPerp2D - rVec%y
                 else
                    lgExit = .true.
                    return
                 end if
              end if
           case ("X")
              xP = 1
              lgHasInPlaneVel = (abs(vHat%y) > 1.e-20 .or. abs(vHat%z) > 1.e-20)
              if (rVec%x < minPerp2D) then
                 if ((lgSymmetricXYZ .or. lgSpecular2D) .and. lgHasInPlaneVel) then
                    vHat%x = -vHat%x
                    rVec%x = 2.0*minPerp2D - rVec%x
                 else
                    lgExit = .true.
                    return
                 end if
              else if (rVec%x > maxPerp2D) then
                 if (lgSpecular2D .and. lgHasInPlaneVel) then
                    vHat%x = -vHat%x
                    rVec%x = 2.0*maxPerp2D - rVec%x
                 else
                    lgExit = .true.
                    return
                 end if
              end if
           end select
        end if

        ! x-axis
        if (.not. (lg2D .and. perpAxis2D == "X")) then
           if (xP < grid(1)%nx) then
              if (rVec%x > (grid(1)%xAxis(xP) + grid(1)%xAxis(xP+1))/2.) then
                 xP = xP + 1
              end if
           else if (xP == grid(1)%nx) then
              if (rVec%x > grid(1)%xAxis(xP)) then
                 lgExit = .true.
                 return
              end if
           end if
           if (xP > 1) then
              if (rVec%x < (grid(1)%xAxis(xP) + grid(1)%xAxis(xP-1))/2.) then
                 xP = xP - 1
              end if
           end if
        end if

        ! y-axis
        if (.not. (lg2D .and. perpAxis2D == "Y")) then
           if (yP < grid(1)%ny) then
              if (rVec%y > (grid(1)%yAxis(yP) + grid(1)%yAxis(yP+1))/2.) then
                 yP = yP + 1
              end if
           else if (yP == grid(1)%ny) then
              if (rVec%y > grid(1)%yAxis(yP)) then
                 lgExit = .true.
                 return
              end if
           end if
           if (yP > 1) then
              if (rVec%y < (grid(1)%yAxis(yP) + grid(1)%yAxis(yP-1))/2.) then
                 yP = yP - 1
              end if
           end if
        end if

        ! z-axis
        if (.not. (lg2D .and. perpAxis2D == "Z")) then
           if (zP < grid(1)%nz) then
              if (rVec%z > (grid(1)%zAxis(zP) + grid(1)%zAxis(zP+1))/2.) then
                 zP = zP + 1
              end if
           else if (zP == grid(1)%nz) then
              if (rVec%z > grid(1)%zAxis(zP)) then
                 lgExit = .true.
                 return
              end if
           end if
           if (zP > 1) then
              if (rVec%z < (grid(1)%zAxis(zP) + grid(1)%zAxis(zP-1))/2.) then
                 zP = zP - 1
              end if
           end if
        end if

        if (.not. lgSymmetricXYZ) then
           if ((xP < 1) .or. (zP < 1) .or. (yP < 1)) then
              lgExit = .true.
              return
           end if
        end if

        ! check if there's any symmetries and if the path must be reflected
        if (lgSymmetricXYZ) then
           if (.not. (lg2D .and. perpAxis2D == "X")) then
              if (xP < 1) then
                 vHat%x = -vHat%x
                 rVec%x = -rVec%x
                 call locate(grid(1)%xAxis, rVec%x, xP)
                 if (xP < grid(1)%nx) then
                    if (rVec%x > (grid(1)%xAxis(xP) + grid(1)%xAxis(xP+1))/2.) xP = xP + 1
                 else if (xP == grid(1)%nx) then
                    if (rVec%x > grid(1)%xAxis(xP)) then
                       lgExit = .true.
                       return
                    end if
                 end if
              end if
           end if
           if (.not. (lg2D .and. perpAxis2D == "Y")) then
              if (yP < 1) then
                 vHat%y = -vHat%y
                 rVec%y = -rVec%y
                 call locate(grid(1)%yAxis, rVec%y, yP)
                 if (yP < grid(1)%ny) then
                    if (rVec%y > (grid(1)%yAxis(yP) + grid(1)%yAxis(yP+1))/2.) yP = yP + 1
                 else if (yP == grid(1)%ny) then
                    if (rVec%y > grid(1)%yAxis(yP)) then
                       lgExit = .true.
                       return
                    end if
                 end if
              end if
           end if
           if (.not. (lg2D .and. perpAxis2D == "Z")) then
              if (zP < 1) then
                 vHat%z = -vHat%z
                 rVec%z = -rVec%z
                 call locate(grid(1)%zAxis, rVec%z, zP)
                 if (zP < grid(1)%nz) then
                    if (rVec%z > (grid(1)%zAxis(zP) + grid(1)%zAxis(zP+1))/2.) zP = zP + 1
                 else if (zP == grid(1)%nz) then
                    if (rVec%z > grid(1)%zAxis(zP)) then
                       lgExit = .true.
                       return
                    end if
                 end if
              end if
           end if
        end if

    end subroutine stepPathTauCell


    ! This subroutine performs the integration of opacities from a position
    ! aVec to the edge of the nebula in the arbitrary direction uHat
    subroutine integratePathTau(grid, aVec, uHat, freq, nTau, absTau, lambda)
        implicit none

        type(grid_type), intent(in) :: grid(*)        ! grid

        type(vector),  intent(in)   :: uHat           ! direction vector
        type(vector),  intent(in)   :: aVec           ! starting position vector

        integer, intent(out)        :: nTau           ! actual size of opacity and length arrays

        real, intent(in)            :: freq           ! the frequency [Hz]
        real, intent(out), &
             & dimension(maxTau)    :: absTau         ! absorption optical depth

        real, intent(out), &
             & dimension(maxTau)    :: lambda         ! distance array

        ! local variables

        type(vector)                :: rVec           ! position vector
        type(vector)                :: vHat           ! direction vector

        integer                     :: freqP          ! frequency index
        integer                     :: i              ! counter
        integer                     :: xP, yP, zP     ! x, y, and z axis indeces
        integer, parameter          :: nl = 4         ! defines the size of the distance increment

        real                        :: dlSmall        ! distance increment
        logical                     :: lgExit

        if (nGrids>1) then
           print*, 'integratePathTau: tau integration routine still not developed for multiple grids... returning '
           return
        end if

        ! locate the input frequency in the frequency array
        call locate(nuArray(1:nbins), freq, freqP)
        if (freqP<=0 .or. freqP>nbins) then
           print*, "! integratePathTau: out of bounds frequency ", freqP, freq, nuArray(1), nuArray(nbins)
           print*, nuArray
           stop
        end if

        call initPathTauStart(grid, aVec, xP, yP, zP, "! integratePathTau")

        ! define dlSmall
        dlSmall = dl(1)

        ! initialize the optical depth arrays
        absTau = 0.
        lambda = 0.

        ! initialize nTau
        nTau = 1

        ! initialize direction vactor
        vHat = uHat

        ! the first optical depth are all zero (as displacement is zero)
        ! these have already been set to zero in the initialization of the arrays
        ! so add the first displacements vector (dlSmall*vHat)
        rVec = aVec + dlSmall*vHat

        ! execute this loop until the edge of the grid is reached
        do i = 1, maxTau

           call stepPathTauCell(grid, rVec, vHat, xP, yP, zP, lgExit)
           if (lgExit) exit

           ! increment nTau
           nTau = nTau + 1

           ! this should never happen (as the loop is stopped at maxTau anyway)
           ! just a sanity check
           if (nTau > maxTau) then
               print*, "! integratePathTau: tau arrays are full", xP,yP,zP
               exit
           end if

           ! add the delta optical depths
           absTau(nTau) = absTau(nTau-1) + grid(1)%opacity(grid(1)%active(xP,yP,zP),freqP)*dlSmall
           lambda(nTau) = lambda(nTau-1) + dlSmall
           ! increment position by adding the displacements vector (dlSmall*vHat)
           rVec = rVec + dlSmall*vHat

        end do

    end subroutine integratePathTau


    ! This subroutine performs the integration of opacities from a position
    ! aVec to the edge of the nebula in the arbitrary direction uHat
    subroutine integratePathTauNu(grid, aVec, uHat, freqP, absTau)
        implicit none

        type(grid_type), intent(in) :: grid(*)        ! grid

        type(vector),  intent(in)   :: uHat           ! direction vector
        type(vector),  intent(in)   :: aVec           ! starting position vector

        real, intent(out)           :: absTau         ! extinction optical depth

        ! local variables

        type(vector)                :: rVec           ! position vector
        type(vector)                :: vHat           ! direction vector

        integer, intent(in)         :: freqP           ! the frequency index

        integer                     :: i              ! counter
        integer                     :: xP, yP, zP     ! x, y, and z axis indeces

        real                        :: dlSmall        ! distance increment
        logical                     :: lgExit

        if (nGrids>1) then
           print*, 'integratePathTau: tau integration routine still not developed for multiple grids... returning '
           return
        end if

        call initPathTauStart(grid, aVec, xP, yP, zP, "! integratePathTauNu")

        dlSmall = 0.
        if (lg2D .and. lgCylinder2D) then
           ! Radial step along xAxis applies to both x and y velocity components in cylindrical (R, Z)
           if ((abs(uHat%x) > 0. .or. abs(uHat%y) > 0.) .and. grid(1)%nx > 1) then
              dlSmall = abs(grid(1)%xAxis(2) - grid(1)%xAxis(1))
              do i = 2, grid(1)%nx-1
                 dlSmall = min(dlSmall, abs(grid(1)%xAxis(i+1) - grid(1)%xAxis(i)))
              end do
           end if
           if (abs(uHat%z) > 0. .and. grid(1)%nz > 1) then
              if (dlSmall <= 0.) then
                 dlSmall = abs(grid(1)%zAxis(2) - grid(1)%zAxis(1))
                 do i = 2, grid(1)%nz-1
                    dlSmall = min(dlSmall, abs(grid(1)%zAxis(i+1) - grid(1)%zAxis(i)))
                 end do
              else
                 do i = 1, grid(1)%nz-1
                    dlSmall = min(dlSmall, abs(grid(1)%zAxis(i+1) - grid(1)%zAxis(i)))
                 end do
              end if
           end if
        else
           if (uHat%x > 0. .and. grid(1)%nx > 1) then
              ! find linear increment
              dlSmall =  abs(grid(1)%xAxis(2) - grid(1)%xAxis(1))
              do i = 2, grid(1)%nx-1
                 dlSmall = min(dlSmall, abs(grid(1)%xAxis(i+1)-grid(1)%xAxis(i)) )
              end do
           end if

           if (uHat%y > 0. .and. grid(1)%ny > 1) then
              if (dlSmall <= 0.) then
                 ! find linear increment
                 dlSmall =  abs(grid(1)%yAxis(2) - grid(1)%yAxis(1))
                 do i = 2, grid(1)%ny-1
                    dlSmall = min(dlSmall, abs(grid(1)%yAxis(i+1)-grid(1)%yAxis(i)) )
                 end do
              else
                 ! find linear increment
                 do i = 1, grid(1)%ny-1
                    dlSmall = min(dlSmall, abs(grid(1)%yAxis(i+1)-grid(1)%yAxis(i)) )
                 end do
              end if
           end if

           if (uHat%z > 0. .and. grid(1)%nz > 1) then
              if (dlSmall <= 0.) then
                 ! find linear increment
                 dlSmall =  abs(grid(1)%zAxis(2) - grid(1)%zAxis(1))
                 do i = 2, grid(1)%nz-1
                    dlSmall = min(dlSmall, abs(grid(1)%zAxis(i+1)-grid(1)%zAxis(i)) )
                 end do
              else
                 ! find linear increment
                 do i = 1, grid(1)%nz-1
                    dlSmall = min(dlSmall, abs(grid(1)%zAxis(i+1)-grid(1)%zAxis(i)) )
                 end do
              end if
           end if
        end if

        if (dlSmall <= 0.) dlSmall = dl(1)
        dlSmall = dlSmall/2.

        ! initialize the optical depth arrays
        absTau = 0.

        ! initialize direction vector
        vHat = uHat

        ! the first optical depth are all zero (as displacement is zero)
        ! these have already been set to zero in the initialization of the arrays
        ! so add the first displacements vector (dlSmall*vHat)
        rVec = aVec + dlSmall*vHat

        ! execute this loop until the edge of the grid is reached
        do i = 1, maxTau

           call stepPathTauCell(grid, rVec, vHat, xP, yP, zP, lgExit)
           if (lgExit) exit

           ! add the delta optical depths
           absTau = absTau + grid(1)%opacity(grid(1)%active(xP,yP,zP),freqP)*dlSmall

           ! increment position by adding the displacements vector (dlSmall*vHat)
           rVec = rVec + dlSmall*vHat

        end do

    end subroutine integratePathTauNu

end module pathIntegration_mod
