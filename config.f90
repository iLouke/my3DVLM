module config
    implicit none
    type Configuration
        ! Parameters
        real(8) :: Vinit        = 0.d0
        real(8) :: yaw_deg      = 0.d0
        real(8) :: pitch_deg    = 0.d0
        real(8) :: roll_deg     = 0.d0
        real(8) :: epsilon      = 0.d0
        real(8) :: flheight     = 0.d0
        real(8) :: dx           = 0.d0 ! Biggest edge of all panels. This is used to calculate dt !!!EXPERIMENTAL!!!
        real(8) :: dt           = 0.d0
        ! Transformed
        real(8) :: yaw_rad      = 0.d0
        real(8) :: pitch_rad    = 0.d0
        real(8) :: roll_rad     = 0.d0
        ! Indexes
        integer :: symmetry     = 0
        integer :: inches       = 0
        integer :: grnd_eff     = 0
        ! Hardcode the filename for better connection with pre/post processors
        character(10) :: panelsf    = 'panels.dat'
        character(10) :: pointsf    = 'points.dat'
        ! Number of points and panels for allocations
        integer :: npoints = 0, npanels = 0, nmark = 0

        contains

        procedure :: set_angles_to_rad   => set_angles_to_rad
        procedure :: set_inches_to_meter => set_inches_to_meter
    end type Configuration

    contains

    subroutine set_inches_to_meter(this)
        implicit none
        class(Configuration) :: this
        this%flheight = inch2meter(this%flheight)
    end subroutine set_inches_to_meter

    subroutine set_angles_to_rad(this)
        implicit none
        class(Configuration) :: this
        this%yaw_rad   = deg2rad(this%yaw_deg)
        this%pitch_rad = deg2rad(this%pitch_deg)
        this%roll_rad  = deg2rad(this%roll_deg)
    end subroutine set_angles_to_rad

    function deg2rad(angle_deg) result(angle_rad)
        implicit none
        real(8) :: angle_deg
        real(8) :: angle_rad
        real(8), parameter :: pi = 4.d0*atan(1.d0)
        angle_rad  = angle_deg * pi / 180.d0
    end function

    function inch2meter(length_in) result(length_met)
        implicit none
        real(8) :: length_in
        real(8) :: length_met
        length_met = length_in * 0.0254
    end function inch2meter

end module config