program main
    use panel_3d
    use config
    implicit none
    ! Constants
    real(8), parameter :: pi  = 4.d0*atan(1.d0)
    real(8), parameter :: rho = 1.225
    ! Types
    type(Panel_), allocatable :: panel(:)
    type(Configuration)       :: setting
    ! Allocatable Arrays
    real(8), allocatable :: A  (:,:)    ! AIC Matrix            (N_panels **2)
    real(8), allocatable :: DW (:,:)    ! Downwash Matrix       ()

    ! Variables for Computations
    real(8) :: Lift = 0.d0, Drag = 0.d0
    real(8) :: MomentX = 0.d0, MomentY = 0.d0, MomentZ = 0.d0
    real(8) :: CL = 0.d0, CD = 0.d0
    real(8) :: CMx = 0.d0, CMy = 0.d0, CMz = 0.d0
    real(8) :: time = 0.d0
    ! Counters
    integer :: iteration = 0, max_iterations = 100   
    integer :: i,j,k,l,m
    ! Variables for System
    real    :: CPU_start = 0., CPU_end = 0., CPU_run = 0.

    ! Gets the start time of solver
    call cpu_time(CPU_start)
    ! Program Starts Here

    ! Reading of inputs
    call read_inputs(setting)
    ! As the input file is read, we allocate memory needed for panels
    allocate(panel(setting%npanels))
    ! Reading of geometry
    call read_geometry(setting,panel)
    ! Set timestep (dt) as dt = dx / vinit / 4.d0 
    ! where dx is the biggest panel edge found in geometry mesh (calculated in read_geometry)
    call timestep_calc()

















    call cpu_time(CPU_end)
    CPU_run = CPU_end - CPU_run
    write (6,100) CPU_run/60., CPU_run  ! Writes on screen the runtime of solver

    ! FORMATS
    100 format('Solver was executed in : ',f3.1,' m ( ',f5.2,' s )')

    contains

    subroutine timestep_calc()
        setting%dt = setting%dx / setting%Vinit / 4.d0
        time = - setting%dt
    end subroutine timestep_calc
end program main

subroutine read_inputs(set)
    use config
    implicit none
    type(Configuration), intent(inout) :: set
    character(65) :: fmt
    integer :: count_lines

    open(unit = 10, file = 'inputs.dat', status = 'unknown')

    fmt = '(F10.3,/,F10.3,/,F10.3,/,F10.3,/,F10.3,/,F10.3,/,I10,/,I10,/,I10)'

    read(10,fmt) set%Vinit, set%pitch_deg, set%yaw_deg, set%roll_deg, &
    &            set%epsilon, set%flheight, set%symmetry, set%inches, &
    &            set%grnd_eff
    
    close(10)

    call set%set_angles_to_rad()
    ! If the inches = 1 then switch them to meters
    if (set%inches == 1) call set%set_inches_to_meter()
    ! Count the points and panels
    set%npoints = count_lines(set%pointsf)
    set%npanels = count_lines(set%panelsf)
end subroutine read_inputs

subroutine read_geometry(set,p)
    use panel_3d
    use config
    implicit none
    type(Configuration), intent(inout)     :: set
    type(Panel_)       , intent(inout)  :: p(set%npanels)
    ! Local Variables
    integer :: i, j, npan, npoi
    real(8) , allocatable :: x(:), y(:), z(:)
    real(8)               :: dist = 0.d0, max_dist = 0.d0
    integer , allocatable :: mark(:)
    integer :: p1, p2 ,p3 ,p4, sec
    integer :: c1,c2,c

    npoi = set%npoints
    npan = set%npanels

    ! Allocate temp arrays
    allocate(x(npoi))
    allocate(y(npoi))
    allocate(z(npoi))
    allocate(mark(npoi))

    ! Reading Points
    open(unit = 10, file = set%pointsf, status = 'unknown')

    do i = 1, npoi
        read(10,100) j, x(i), y(i), z(i), mark(i) 

        ! Add flight height to z coordinate
        if (set%flheight > 0.d0) z(i) = z(i) + set%flheight
    end do
    ! Update Number of marked vertices
    set%nmark = sum(mark)
    close(10)
    
    ! Reading Panels
    open(unit = 20, file = set%panelsf, status = 'unknown')

    do i = 1, npan
        read(20,200) j, p1, p2, p3, p4, sec

        ! Saving point ID
        p(i)%pid = [p1,p2,p3,p4]
        ! Check if its triangle
        if (p3==p4) p(i)%is_triangle = .true.
        ! Saving Corresponding points
        p(i)%x = [x(p1),x(p2),x(p3),x(p4)]
        p(i)%y = [y(p1),y(p2),y(p3),y(p4)]
        p(i)%z = [z(p1),z(p2),z(p3),z(p4)]
        ! Saving Section that the panel belongs to
        p(i)%section = sec
        ! Save if a point of the panel is marked
        if (sum([mark(p1),mark(p2),mark(p3),mark(p4)]) > 0) then
            p(i)%mark = [mark(p1),mark(p2),mark(p3),mark(p4)]
        end if
        ! ====================================
        ! Calculate the geometrical properties
        ! ====================================
        if (set%inches == 1) call p(i)%convert_to_meters()
        ! Rotate panel(i)%x ,y,z by yaw,pitch,roll
        call p(i)%rotate_panel(set%pitch_rad,set%yaw_rad,set%roll_rad)
        ! Sets panel(i)%cp 
        call p(i)%get_collocation_point()
        ! Sets panel(i)%normal
        call p(i)%get_normal_vector()
        ! Sets panel(i)%ds
        call p(i)%get_panel_surface()
        ! ====================================
        ! Calculate biggest edge of each panel
        ! ====================================
        do c = 1, 4
            c1 = c
            c2 = c + 1
            if (c2 == 5) c2 = 1

            dist = sqrt( (x(c1)-x(c2))**2 + &
            &            (y(c1)-y(c2))**2 + &
            &            (z(c1)-z(c2))**2 )

            if (max_dist < dist) max_dist = dist
        end do
    end do
    ! Set dx as the biggest edge of panel (dx is used to calculate the timestep)
    set%dx = max_dist

    close(20)

    ! Deallocate temp arrays just to be sure
    deallocate(x,y,z)
    deallocate(mark)

    ! Formats
    100 format (I10,3f15.5,I10)
    200 format (6I10)
end subroutine

function count_lines(filename) result(nlines)
    implicit none
    character(len=*)    :: filename
    integer             :: nlines
    integer             :: io
    open(10,file=filename, iostat=io, status='old')
    if (io/=0) stop 'Cannot open file! '

    nlines = 0
    do
        read(10,*,iostat=io)
        if (io/=0) exit
        nlines = nlines + 1
    end do

    close(10)
end function count_lines