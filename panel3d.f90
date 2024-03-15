module panel_3d
    implicit none
    type Panel_
        integer :: pid(4)       = 0         ! Points ID (from points file) consisting Panel (May be useless)
        real(8) :: x(4)         = 0.d0      ! The xyz of the four points consisting a quad panel.
        real(8) :: y(4)         = 0.d0      ! If point 3 is equal to point 4 then its a triangle
        real(8) :: z(4)         = 0.d0      ! ==================================================
        real(8) :: cp(3)        = 0.d0      ! Collocation Point of the panel (centroid)
        real(8) :: norm(3)      = 0.d0      ! Normal Vector of the panel
        real(8) :: ds           = 0.d0      ! Surface of the panel
        logical :: is_triangle  = .false.   ! Boolean value if the panel is triangle
        integer :: mark(4)      = 0
        integer :: section      = 0

        contains
        procedure :: convert_to_meters      => convert_to_meters
        procedure :: rotate_panel           => rotate_panel
        procedure :: get_collocation_point  => get_collocation_point
        procedure :: get_normal_vector      => get_normal_vector
        procedure :: get_panel_surface      => get_panel_surface
    end type Panel_

    contains

    subroutine convert_to_meters(this)
        implicit none
        class(Panel_) :: this
        integer :: i
        
        do i = 1, 4
            this%x(i) = inch2meter(this%x(i))
            this%y(i) = inch2meter(this%y(i))
            this%z(i) = inch2meter(this%z(i))
        end do

    end subroutine convert_to_meters

    subroutine rotate_panel(this,alpha,beta,gamma)
        implicit none
        class(Panel_) :: this
        real(8), intent(in):: alpha,beta,gamma
        real(8) :: sina,cosa,sinb,cosb,sing,cosg
        real(8) :: tempx(4),tempy(4),tempz(4)

        sina = sin(alpha); cosa = cos(alpha)
        sinb = sin(beta) ; cosb = cos(beta)
        sing = sin(gamma); cosg = cos(gamma)

        tempx = this%x*cosa*cosb + this%y*sinb + this%z*sina*cosb

        tempy = - this%x*(cosa*sinb*cosg - sina*sing) &
                + this%y*(cosb*cosg                 ) &
                - this%z*(sina*sinb*cosg - cosa*sing)

        tempz =   this%x*(cosa*sinb*sing - sina*cosg) &
        &       - this%y*(cosb*sing                 ) &
        &       + this%z*(sina*sinb*sing + cosa*cosg)

        this%x = tempx
        this%y = tempy
        this%z = tempz
    end subroutine rotate_panel

    subroutine get_collocation_point(this)
        implicit none
        class(Panel_) :: this

        if (this%is_triangle) then
            this%cp(1) = sum(this%x(1:3))/3.d0
            this%cp(2) = sum(this%y(1:3))/3.d0
            this%cp(3) = sum(this%z(1:3))/3.d0
        else
            this%cp(1) = sum(this%x(1:4))/4.d0
            this%cp(2) = sum(this%y(1:4))/4.d0
            this%cp(3) = sum(this%z(1:4))/4.d0
        end if
    end subroutine get_collocation_point

    subroutine get_normal_vector(this)
        implicit none
        class(Panel_) :: this
        real(8) :: p1(3),p2(3),p3(3)

        p1 = [this%x(1),this%y(1),this%z(1)]
        p2 = [this%x(2),this%y(2),this%z(2)]
        p3 = [this%x(3),this%y(3),this%z(3)]

        call plane_exp_normal_3d (p1,p2,p3,this%norm)
    end subroutine get_normal_vector

    subroutine get_panel_surface(this)
        implicit none
        class(Panel_) :: this
        real(8) :: q(3,4)

        if (this%is_triangle) then
            ! Sets the input to the libgeom.f90 subroutine
            q(1,1) = this%x(1); q(1,2) = this%x(2); q(1,3) = this%x(3)
            q(2,1) = this%y(1); q(2,2) = this%y(2); q(2,3) = this%y(3)
            q(3,1) = this%z(1); q(3,2) = this%z(2); q(3,3) = this%z(3)

            call triangle_area_3d(q,this%ds)
        else
            ! Sets the input to the libgeom.f90 subroutine
            q(1,1) = this%x(1); q(1,2) = this%x(2); q(1,3) = this%x(3); q(1,4) = this%x(4);
            q(2,1) = this%y(1); q(2,2) = this%y(2); q(2,3) = this%y(3); q(2,4) = this%y(4);
            q(3,1) = this%z(1); q(3,2) = this%z(2); q(3,3) = this%z(3); q(3,4) = this%z(4);

            call quad_area_3d(q,this%ds)
        end if
    end subroutine get_panel_surface

    function inch2meter(length_in) result(length_met)
        implicit none
        real(8) :: length_in
        real(8) :: length_met
        length_met = length_in * 0.0254
    end function inch2meter

end module panel_3d