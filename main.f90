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
    real(8), allocatable :: A     (:,:)    ! AIC Matrix            (N_panels **2)
    real(8), allocatable :: B     (:,:)    ! DownWash              (N_panels **2)
    real(8), allocatable :: GAMA  (:)      ! Vortex Strength       (N_panels)
    real(8), allocatable :: RHS   (:)      ! RHS                   (N_panels)
    real(8), allocatable :: GAMAw (:,:)    ! Wake Panel Strength   (max_iter, N_panels_wake)

    ! Variables for Computations
    real(8) :: Lift = 0.d0, Drag = 0.d0
    real(8) :: MomentX = 0.d0, MomentY = 0.d0, MomentZ = 0.d0
    real(8) :: CL = 0.d0, CD = 0.d0
    real(8) :: CMx = 0.d0, CMy = 0.d0, CMz = 0.d0
    real(8) :: time = 0.d0
    ! Counters
    integer :: iteration = 0, max_iterations = 100   
    integer :: i,j,k,l,m, cmax, c, p1,p2
    ! Variables for System
    real    :: CPU_start = 0., CPU_end = 0., CPU_run = 0.

    ! Gets the start time of solver
    call cpu_time(CPU_start)
    ! Program Starts Here

    ! Reading of inputs
    call read_inputs(setting)
    ! As the input file is read, we allocate memory needed for panels
    allocate(panel(setting%npanels))
    allocate(A(setting%npanels, setting%npanels))
    allocate(B(setting%npanels, setting%npanels))
    allocate(GAMA(setting%npanels))
    allocate(RHS(setting%npanels))
    ! Reading of geometry
    call read_geometry(setting,panel)
    ! Set timestep (dt) as dt = dx / vinit / 4.d0 
    ! where dx is the biggest panel edge found in geometry mesh (calculated in read_geometry)
    call timestep_calc()

    ! Starting Loop
    do iteration = 1, max_iterations

        call vorcalc(iteration,setting%npanels,panel, &
        & A,B,GAMA,RHS,setting%symmetry,setting%grnd_eff,setting%vinit)


    end do
    
    
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

subroutine vorcalc(iter,n,p,A,B,GAMA,RHS,sym,grnd,vinit)
    use panel_3d
    implicit none
    integer     , intent(in)   :: iter, n, sym, grnd
    real(8)     , intent(in)   :: vinit
    type(Panel_), intent(in)   :: p(n)
    real(8)     ,intent(inout) :: A(n,n),B(n,n),GAMA(n),RHS(n)
    integer :: i,j,c,c1,c2, cmax
    real(8) ::       u,      v,      w
    real(8) ::    sumU,   sumV,   sumW
    real(8) ::  s_sumU, s_sumV, s_sumW
    real(8)  :: w_sumU, w_sumV, w_sumW
    real(8) :: x,y,z, x1,x2,y1,y2,z1,z2

    GAMA = 0.d0; RHS = 0.d0

    ! If it is the first iteration then the A matrix should be computed
    if (iter == 1) then
        A = 0.d0; B = 0.d0
        do i = 1,n
            x = p(i)%cp(1)
            y = p(i)%cp(2)
            z = p(i)%cp(3)
            do j = 1,n

                  sumU = 0.d0;   sumV = 0.d0;   sumW = 0.d0
                s_sumU = 0.d0; s_sumV = 0.d0; s_sumW = 0.d0
                if (p(j)%is_triangle) then
                    cmax = 3
                else
                    cmax = 4
                end if


                do c = 1, cmax
                    c1 = c
                    c2 = c + 1
                    if (c2 == cmax + 1) c2 =1

                    x1 = p(j)%x(c1); y1 = p(j)%y(c1); z1 = p(j)%z(c1);
                    x2 = p(j)%x(c2); y2 = p(j)%y(c2); z2 = p(j)%z(c2);

                    call vortex(x,y,z,x1,y2,z1,x2,y2,z2,1.d0,u,v,w)

                    sumU = sumU + u; sumV = sumV + v; sumW = sumW + w;
                end do
                ! Symmetry Half
                if (sym==1) then
                    do c = 1, cmax
                        c1 = c
                        c2 = c + 1
                        if (c2 == cmax + 1) c2 =1
                        x1 = p(j)%x(c1); y1 = p(j)%y(c1); z1 = p(j)%z(c1);
                        x2 = p(j)%x(c2); y2 = p(j)%y(c2); z2 = p(j)%z(c2);
    
                        call vortex(x,-y,z,x1,y2,z1,x2,y2,z2,1.d0,u,v,w)
    
                        sumU = sumU + u; sumV = sumV - v; sumW = sumW + w;
                    end do
                end if
                ! Ground Mirror Image
                if (grnd==1) then
                    do c = 1, cmax
                        c1 = c
                        c2 = c + 1
                        if (c2 == cmax + 1) c2 =1
                        x1 = p(j)%x(c1); y1 = p(j)%y(c1); z1 = p(j)%z(c1);
                        x2 = p(j)%x(c2); y2 = p(j)%y(c2); z2 = p(j)%z(c2);
    
                        call vortex(x,y,-z,x1,y2,z1,x2,y2,z2,1.d0,u,v,w)
    
                        sumU = sumU + u; sumV = sumV + v; sumW = sumW - w;
                    end do
                end if
                ! Symmetry Ground Mirror Image
                if (grnd==1) then
                    do c = 1, cmax
                        c1 = c
                        c2 = c + 1
                        if (c2 == cmax + 1) c2 =1
                        x1 = p(j)%x(c1); y1 = p(j)%y(c1); z1 = p(j)%z(c1);
                        x2 = p(j)%x(c2); y2 = p(j)%y(c2); z2 = p(j)%z(c2);
    
                        call vortex(x,-y,-z,x1,y2,z1,x2,y2,z2,1.d0,u,v,w)
    
                        sumU = sumU + u; sumV = sumV - v; sumW = sumW - w;
                    end do
                end if

                ! Sum up all the induced velocities
                A(i,j) = A(i,j) + p(i)%norm(1)*sumU + p(i)%norm(2)*sumV + p(i)%norm(3)*sumW
            end do
            RHS(i) = - (Vinit*p(i)%norm(1))
        end do
    else
        do i = 1, n
            x = p(i)%cp(1)
            y = p(i)%cp(2)
            z = p(i)%cp(3)

            w_sumU = 0.d0; w_sumV = 0.d0; w_sumW = 0.d0
            
            ! Add wake influence
            print *, 'NO WAKE INFLUENCE CALC IN VORCALC'


            RHS(i) = -(    Vinit*p(i)%norm(1) &
            &           + w_sumU*p(i)%norm(1) &
            &           + w_sumV*p(i)%norm(2) &
            &           + w_sumW*p(i)%norm(3) )
        end do
    end if

    ! Solve for Gamma
    GAMA = RHS
    call solve_axb(A,n,GAMA)
end subroutine vorcalc

subroutine solve_axb(a, n, b)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n, n), intent(in)    :: a
    real(8), dimension(n)   , intent(inout) :: b
  
    integer :: info
    real(8), dimension(n) :: ipiv(n)
    ! Call DGESV to solve the system AX = B. X is saved in B
    call dgesv(n,1,a,n,ipiv,b,n,info)
    if (info .ne. 0) then
      print *, 'Error in dgesv: ', info
      ! Handle the error appropriately
    end if
end subroutine solve_axb


subroutine vortex(X,Y,Z,X1,Y1,Z1,X2,Y2,Z2,GAMA,U,V,W)
    implicit double precision (A-Z)
    ! SUBROUTINE VORTEX CALCULATES THE INDUCED VELOCITY (U,V,W) AT A POI
    ! (X,Y,Z) DUE TO A VORTEX ELEMENT VITH STRENGTH GAMA PER UNIT LENGTH
    ! POINTING TO THE DIRECTION (X2,Y2,Z2)-(sX1,Y1,Z1).
    PAY = 4.d0*atan(1.d0)
    RCUT=1.0E-10
    ! CALCULATION OF R1 X R2
    R1R2X=(Y-Y1)*(Z-Z2)-(Z-Z1)*(Y-Y2)
    R1R2Y=-((X-X1)*(Z-Z2)-(Z-Z1)*(X-X2))
    R1R2Z=(X-X1)*(Y-Y2)-(Y-Y1)*(X-X2)
    ! CALCULATION OF (R1 X R2 )**2
    SQUARE=R1R2X*R1R2X+R1R2Y*R1R2Y+R1R2Z*R1R2Z
    ! CALCULATION OF R0(R1/R(R1)-R2/R(R2))
    R1=SQRT((X-X1)*(X-X1)+(Y-Y1)*(Y-Y1)+(Z-Z1)*(Z-Z1))
    R2=SQRT((X-X2)*(X-X2)+(Y-Y2)*(Y-Y2)+(Z-Z2)*(Z-Z2))
    IF((R1.LT.RCUT).OR.(R2.LT.RCUT).OR.(SQUARE.LT.RCUT)) GOTO 1 !GROUND
    R0R1=(X2-X1)*(X-X1)+(Y2-Y1)*(Y-Y1)+(Z2-Z1)*(Z-Z1)
    R0R2=(X2-X1)*(X-X2)+(Y2-Y1)*(Y-Y2)+(Z2-Z1)*(Z-Z2)
    COEF=GAMA/(4.0*PAY*SQUARE)*(R0R1/R1-R0R2/R2)
    U=R1R2X*COEF
    V=R1R2Y*COEF
    W=R1R2Z*COEF
    GOTO 2
    ! WHEN POINT (X,Y,Z) LIES ON VORTEX ELEMENT; ITS INDUCED VELOCITY IS
    1 U=0.
    V=0.
    W=0.
    2 CONTINUE
    RETURN
end subroutine vortex