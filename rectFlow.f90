! ==========================================
!  Program:     /Users/mandli/Documents/School/grad/ns_course
!  File:        projection_method
!  Created:     2010-05-02
!  Author:      Kyle Mandli
! =========================================
!      Copyright (C) 2010-05-02 Kyle Mandli <mandli@amath.washington.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD)
!  license
!                     http://www.opensource.org/licenses/
 ! ======================================= =

module grid_module

    implicit none

    ! Grid parameters
    integer :: N_y , N_x, nlx, nly, n_lagrangian_points
    double precision :: L_x, h, L_y, HL_y, HL_x, H_xOffset, h_l
    character*5 :: fileSuffix

    ! Grid arrays
    ! i-1,j  |  i,j  |  i+1,j           u <-- x_edge, matches p indices
    ! ---v---|---v---|---v---
    !        |       |
    !  i-1,j |  i,j  |  i+1,j
    !    p   u   p   u   p    <--- x_center, y_center
    !        |       |
    ! i-1,j-1| i,j-1 | i+1,j-1
    ! ---v---|---v---|---v--- <--- y_edge
    !        |       |
    double precision, allocatable :: x_edge(:),x_center(:), x_lag(:)
    double precision, allocatable :: y_edge(:),y_center(:), y_lag(:)

contains

    ! Setup and allocate all grid quantities
    subroutine setup_grid()

        implicit none

        ! Local
        integer :: i

        ! Grid spacing, must be equal in x and y directions
        ! base off y dimension (Length, n points)
        !h = (L_x - 0.d0) / (N_x)
        h = (L_y - 0.d0) / (N_y)
        N_x = (L_x - 0.d0) / h

        ! Grid array allocation
        allocate(x_edge(0:N_x+1),x_center(0:N_x+1))
        allocate(y_edge(0:N_y+1),y_center(0:N_y+1))

        ! Calculate grid arrays
        forall (i=0:N_x+1)
            x_edge(i) = 0.d0 + i * h
            x_center(i) = 0.d0 + (i-0.5d0) * h
        end forall
        forall (i=0:N_y+1)
            y_edge(i) = 0.d0 + i * h
            y_center(i) = 0.d0 + (i-0.5d0) * h
        end forall

    end subroutine setup_grid

    ! Setup and allocate all lagrangian points
    subroutine setup_lagrangian_points()

        implicit none


        !locals
        double precision :: x_front, x_back
        double precision :: y_top, y_bottom
        integer :: i, running_i ! number of points in the x and y directions

        !h_l = h/2 ! spacing between lagrangian points

        ! determine number of points in x and y directions
        !nlx = HL_x / h_l
        !nly = HL_y / h_l

        ! all lagrangian points
        allocate(x_lag(1:(n_lagrangian_points)))
        allocate(y_lag(1:(n_lagrangian_points)))

        x_front = H_xOffset
        x_back = H_xOffset + HL_x
        y_top = 0.5 * L_y + 0.5 * HL_y
        y_bottom = 0.5 * L_y - 0.5 * HL_y

        ! start lagrangian points at top left corner of
        ! rectangle and move clockwise
        ! define top left corner
        x_lag(1) = x_front
        y_lag(1) = y_top
        do i=2,nlx
            x_lag(i) = x_lag(i-1) + h_l
            y_lag(i) = y_lag(1)
        enddo
        do i=nlx+1,nlx+nly
            x_lag(i) = x_lag(nlx)
            y_lag(i) = y_lag(i-1) - h_l
        enddo
        do i = nlx+nly+1,nlx+nly+nlx
            x_lag(i) = x_lag(i-1) - h_l
            y_lag(i) = y_lag(nlx+nly)
        enddo
        do i=nlx+nly+nlx+1,nlx+nly+nlx+nly
            x_lag(i) = x_lag(nlx+nly+nlx)
            y_lag(i) = y_lag(i-1) + h_l
        enddo

    end subroutine setup_lagrangian_points

    subroutine output_lagrangian_points()

        ! local
        integer :: i
        open(unit=35,file='_output/lagrangian_points'//fileSuffix//'.dat',access='sequential',status='unknown')
        do i=1,nlx*2+nly*2
            write(35,*) x_lag(i), y_lag(i)
        enddo
        close(35)

    end subroutine output_lagrangian_points

    subroutine output_grid_centers()

        ! local
        integer :: i
        open(unit=35,file='_output/x_points'//fileSuffix//'.dat',access='sequential',status='unknown')
        do i=1,N_x
            write(35,*) x_edge(i)
        enddo
        close(35)
        open(unit=35,file='_output/y_points'//fileSuffix//'.dat',access='sequential',status='unknown')
        do i=1,N_y
            write(35,*) y_edge(i)
        enddo
        close(35)
    end subroutine output_grid_centers

    ! Grid stretching function
!    double precision pure function F(zeta)
!        implicit none
!        double precision, intent(in) :: zeta
        !F =  tanh(gamma*L_y) / (gamma*L_y) * cosh(gamma*(L_y - zeta))**2
!        F = L_y
!    end function F

    !delta_hxy function equation 2.12 of jfm 2010
    ! computes delta function from a vector, eg (x1,y1) - (x2,y2)
    double precision pure function delta_hxy(x1,y1,x2,y2)
        implicit none
        double precision, intent(in) :: x1, y1, x2, y2
        delta_hxy = delta_h1(x1,x2) * delta_h1(y1,y2)
    end function delta_hxy

    ! delta_h1 function
    double precision pure function delta_h1(n1, n2) ! matches eqn 21 roma-peskin-berger
        implicit none
        double precision, intent(in) :: n1, n2
        ! local
        double precision :: r
        r = (n1 - n2) / h
        delta_h1 = 1.d0 / h * phi_r(r)
    end function delta_h1

    ! phi_r function
    double precision pure function phi_r(r) ! matches eqn 22 roma-peskin-berger
        implicit none
        double precision, intent(in) :: r
        ! local
        double precision :: output, absr
        absr = abs(r)
        if (absr<=0.5) then
            output = 1.d0/3.d0 * (1 + sqrt(1-3*absr**2))
        else if (absr<=1.5) then
            output = 1.d0/6.d0 * (5 - 3*absr - sqrt(1-3*(1-absr)**2))
        else
            output = 0.d0
        end if
        phi_r = output

    end function phi_r

    ! Output grid
    subroutine output_grid(frame,t,u,v,p)

        implicit none

        ! Arguments
        integer, intent(in) :: frame
        double precision, intent(in) :: t
        double precision, dimension(0:N_x+1,0:N_y+1) :: u,v,p

        ! Locals
        integer :: i,j


! Open output file and write file header:

        if (t==0) then
        open(unit=70,file='_output/UVP_'//fileSuffix//'.dat',access='sequential',status='unknown')
        write(70,*) "Grid Size", N_x, N_y
        write(70,*) ' VARIABLES= "x", "y", "u", "v", "p"'
        write(70,100) t,N_X,N_Y
        endif

 100 FORMAT('ZONE T="t = ',e26.16,'"',' F=POINT, I=',I5,' J=', I5)

        if (t>0) then
        write(70,100) t,N_X,N_Y
        endif

        ! Write out data
        do j=1,N_y
        do i=1,N_x

        ! Reset to zero if exponent is too large
        if (abs(u(i,j)) < 1d-99) u(i,j) = 0.d0
        if (abs(v(i,j)) < 1d-99) v(i,j) = 0.d0
        if (abs(p(i,j)) < 1d-99) p(i,j) = 0.d0

        write(70,"(5e26.16)") x_edge(i),y_edge(j),u(i,j),v(i,j),p(i,j)

        enddo
        enddo

    end subroutine output_grid

    ! Output grid
    subroutine output_force_grid(frame,t,u_star,v_star)

        implicit none

        ! Arguments
        integer, intent(in) :: frame
        double precision, intent(in) :: t
        double precision, dimension(0:N_x+1,0:N_y+1) :: u_star, v_star

        ! Locals
        integer :: i,j


! Open output file and write file header:

        if (t==0) then
        open(unit=77,file='_output/force_grid.dat',access='sequential',status='unknown')
        write(77,*) "Grid Size", N_y, N_x
        write(77,*) ' VARIABLES= "x", "y", "u", "v"'
        write(77,100) t,N_X,N_Y
        endif

 100 FORMAT('ZONE T="t = ',e26.16,'"',' F=POINT, I=',I5,' J=', I5)

        if (t>0) then
        write(77,100) t,N_X,N_Y
        endif

        ! Write out data
        do j=1,N_y
        do i=1,N_x

        ! Reset to zero if exponent is too large
        if (abs(u_star(i,j)) < 1d-99) u_star(i,j) = 0.d0
        if (abs(v_star(i,j)) < 1d-99) v_star(i,j) = 0.d0

        write(77,"(5e26.16)") x_edge(i),y_edge(j),u_star(i,j),v_star(i,j)

        enddo
        enddo

    end subroutine output_force_grid



end module grid_module

! =========================================
!  Program:     /Users/mandli/Documents/School/grad/ns_course
!  File:        main
!  Created:     2010-05-02
!  Author:      Kyle Mandli
! ========================================
!      Copyright (C) 2010-05-02 Kyle Mandli <mandli@amath.washington.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD)
!  license
!                     http://www.opensource.org/licenses/
! =========================================

program main

    use grid_module

    implicit none

    ! ====================================
    ! Solver parameters
    integer, parameter :: MAX_ITERATIONS = 100000
    double precision, parameter :: TOLERANCE = 1d-4, CFL = 0.8d0
    logical, parameter :: write_star = .false.
    integer :: n_steps

    ! ===================================
    ! Physics
    double precision :: U_inf = 1.d0
    double precision :: rho,nu,Re !rho = 1.d0, nu=1.d-3

    ! ===================================
    ! Velocity and pressures
    double precision, allocatable :: u(:,:),v(:,:),p(:,:),u_star(:,:),v_star(:,:)
    double precision, allocatable :: u_old(:,:), v_old(:,:)

    ! ===================================
    ! Locals
    character*20 :: arg
    integer :: i,j,n,m,frame,i_R,j_R, k, boxLen
    double precision :: R,t,dt,a, lagFuSum, lagFvSum
    double precision, allocatable :: Flux_ux(:,:),Flux_uy(:,:),Flux_vy(:,:)
    double precision :: uu_x,uv_y,uv_x,vv_y,u_xx,u_yy,v_xx,v_yy
    double precision, allocatable :: Q(:,:),b(:),cp(:),cm(:), ustar_lagF(:), vstar_lagF(:)
    double precision :: xGrid, yGrid, xLag, yLag
    ! ===================================

    ! Get command line arguments
    if (iargc() /= 2) then
        print *,"Wrong number of command line arguments, expected 2."
        print *,"    Lx - Length of domain in x"
        print *,"    boxLen - Length of box"
        print *, "got ", iargc()
        stop
    else
        call getarg(1,arg)
        read (arg,'(I10)') N_y
        call getarg(2,arg)
        read (arg,'(I10)') boxLen
    endif
    write(fileSuffix, "(I3,'_',I1)") N_y, boxLen
    print *, fileSuffix

!!!!!!!!!!!!!!!!!!
!!  Parameters: !!
!!!!!!!!!!!!!!!!!!

    !N_x=10  !Number of grid points in x-direction
    !N_y = 128   !Number of grid points in y-direction
    L_x = 30 !Length of box in x-direction
    L_y = 30  !Length of box in y-direction
    n_steps = MAX_ITERATIONS/50 !Interval that u,v and p are printed to UVP.dat
    ! Setup grid and arrays
    call setup_grid()
    call output_grid_centers()

    !!! Lagrangian Points
    HL_y = 0.1 * L_y  ! Length of rect bluff y-direction
    HL_x = boxLen*HL_y ! Length of rect bluff in x-direction
    H_xOffset = 0.25 * L_x ! how far along x before bluff starts, bluff will always be
                           ! centered in the domain.
    h_l = 0.5*h ! spacing of lagrangian points
    ! determine number of points in x and y directions
    nlx = HL_x / h_l
    nly = HL_y / h_l
    n_lagrangian_points = 2*nlx + 2*nly
    call setup_lagrangian_points()
    call output_lagrangian_points()
    ! ===================================

    allocate(Flux_ux(1:N_x+1,0:N_y+1))
    allocate(Flux_uy(0:N_x,0:N_y))
    allocate(Flux_vy(0:N_x+1,0:N_y))
    allocate(Q(1:N_x,1:N_y))
    allocate(b(1:N_y),cp(1:N_y),cm(1:N_y))

    ! Calculate matrix coefficients for Poisson solve
    a = 1.d0 / h**2
    forall (j=1:N_y)
        b(j) = - (2.d0 / h**2 + y_center(j) / h**2 * (y_edge(j) + y_edge(j-1)))
        cp(j) = (y_center(j) * y_edge(j)) / h**2
        cm(j) = (y_center(j) * y_edge(j-1)) / h**2
    end forall

    ! Velocity and pressure arrays
    allocate(u(0:N_x+1,0:N_y+1),u_star(0:N_x+1,0:N_y+1))
    allocate(v(0:N_x+1,0:N_y+1),v_star(0:N_x+1,0:N_y+1))
    allocate(p(0:N_x+1,0:N_y+1))
    allocate(u_old(1:N_x,1:N_y),v_old(1:N_x,1:N_y))

    ! velocity at lagrangian points
    allocate(ustar_lagF(1:(n_lagrangian_points)))
    allocate(vstar_lagF(1:(n_lagrangian_points)))

    ! Inital conditions
    u = U_inf
    !u = 0.d0
    !! create ramping u at inflow for starting condition, should be zero at boundary
    !do j=0,N_y+1
    !do i=0,H_xOffset/h
    !do j=0,N_y+1
    !    u(i,j) = U_inf - (h*i/H_xOffset)*U_inf
    !enddo
    !enddo
    v = 0.d0
    p = 0.d0

    nu = 1.d-3 !1.d0/Re
    Re = 1.d0/nu
    !Re = 1000.d0
    !nu = 1.d0/Re
    rho = 1.d0
    dt =  h / (Re * U_inf) !* CFL
    !dt = h / (Re * U_inf)
    !dt = h / 1000.
    t = 0.d0
    frame = 0



    ! Output inital condition
    call output_grid(frame,t,u,v,p)
    print "(a,i3,a,i4,a,e16.8)","Writing frame ",frame," during step n=",0," t=",t

    ! Open up file to store residual information in
    open(unit=13, file='_output/residual_'//fileSuffix//'.dat', status="unknown", action="write")

    ! ===================================
    ! Main algorithm loop
    do n=1,MAX_ITERATIONS
        ! Store old step for convergence test
        u_old = u(1:N_x,1:N_y)
        v_old = v(1:N_x,1:N_y)

    ! ===================================
        ! Apply BC  V = x, U = o, P = +
        call bc(u,v,U_inf)

    ! ===================================
        ! Step 1: Update velocity to intermediate step
        ! Calculate fluxes at each boundary
        !
        !  U-Fluxes
        !                                   |   u   |
        !         FUY_i,j                i,j| i,j+1 |i+1,j
        ! -------|---o---|-------    -------v-------v-------
        !        |       |            i-1,j |  i,j  |  i+1,j
        !    FUX_i,j    FUX_i+1,j       u   |   u   |   u
        !        |       |                  |       |
        ! -------|---o---|-------    -------v-------v-------
        !         FUY_i,j-1            i,j-1| i,j-1 | i+1,j-1
        !                                   |   u   |
        !
        !  V-Fluxes
        !                                   |   v   |
        !         FVY_i,j            i-1,j+1| i,j+1 |i,j+1
        ! -------|---o---|-------    -------u---o---u-------
        !        |       |            i-1,j |  i,j  |  i+1,j
        !    FUY_i-1,j  FUY_i,j         v   o   v   o   v
        !        |       |                  |       |
        ! -------|---o---|-------    -------u---o---u-------
        !         FVY_i,j-1            i-1,j| i,j-1 |i,j
        !                                   |   v   |
        forall (i=1:N_x+1,j=0:N_y+1) Flux_ux(i,j) = (u(i-1,j)+u(i,j))**2 / 4.d0
        forall (i=0:N_x,j=0:N_y) Flux_uy(i,j) = (u(i,j)+u(i,j+1)) * (v(i+1,j)+v(i,j)) / 4.d0
        forall (i=0:N_x+1,j=0:N_y) Flux_vy(i,j) = (v(i,j+1) + v(i,j))**2 / 4.d0
        do j=1,N_y
            do i=1,N_x
                ! Advective terms, see Lecture 4 notes, page 12
                uu_x = (Flux_ux(i+1,j) - Flux_ux(i,j)) / h
                uv_y = y_center(j) * (Flux_uy(i,j) - Flux_uy(i,j-1)) / h !

                uv_x = (Flux_uy(i,j) - Flux_uy(i-1,j)) / h
                vv_y = y_edge(j) * (Flux_vy(i,j) - Flux_vy(i,j-1)) / h

                ! Diffusive terms
                u_xx = (u(i+1,j) - 2*u(i,j) + u(i-1,j)) / (h**2)
                u_yy = (y_center(j) / h**2) * (y_edge(j) * (u(i,j+1) - u(i,j)) - y_edge(j-1) * (u(i,j) - u(i,j-1))) !edge leads center by 0.5, verify by looking at poisson solver
                v_xx = (v(i+1,j) - 2.d0*v(i,j) + v(i-1,j)) / (h**2)
                v_yy = (y_edge(j) / h**2) * (y_center(j+1) * (v(i,j+1) - v(i,j)) - y_center(j) * (v(i,j) - v(i,j-1)))

                ! Update to u* and v* values
                u_star(i,j) = u(i,j) + dt * (-1.d0*(uu_x + uv_y) + nu*(u_xx + u_yy))
                v_star(i,j) = v(i,j) + dt * (-1.d0*(uv_x+vv_y) + nu*(v_xx + v_yy))
            enddo
        enddo

        ! Debug, save out u_star,v_star,p
        if (write_star) then
            frame = frame + 1
            print "(a,i3,a,i4,a,e16.8)","Writing frame ",frame," during step n=",n," t=",t
            call output_grid(frame,t,u_star,v_star,p)
        endif
    !===================================
        ! compute u* and v* forcing at lagrangian points
        do k=1,n_lagrangian_points
            ! loop over whole domain (could likely just do a subset to speed up later)
            ! initialize to zero
            ustar_lagF(k) = 0.d0
            vstar_lagF(k) = 0.d0
            do j=1,N_y
                do i=1,N_x
                    ! convenience variables
                    ! WARNING: use edge or centers???!!!!
                    xGrid = x_edge(i)
                    yGrid = y_edge(j)
                    xLag = x_lag(k)
                    yLag = y_lag(k)
                    ustar_lagF(k) = ustar_lagF(k) + (-1.d0) * u_star(i,j) * delta_hxy(xGrid,yGrid,xLag,yLag)*h**2/dt
                    vstar_lagF(k) = vstar_lagF(k) + (-1.d0) * v_star(i,j) * delta_hxy(xGrid,yGrid,xLag,yLag)*h**2/dt
                enddo
            enddo
        enddo

        ! calcluate the force grid (could be combined with the next step, but I wanna look at the force grid)
        do j=1,N_y
            do i=1,N_x
                do k=1,n_lagrangian_points
                    ! sum up forces acting on this point
                    xGrid = x_edge(i)
                    yGrid = y_edge(j)
                    xLag = x_lag(k)
                    yLag = y_lag(k)
                    u_star(i,j) = u_star(i,j) + ustar_lagF(k) * delta_hxy(xGrid,yGrid,xLag,yLag)*h_l**2*dt
                    v_star(i,j) = v_star(i,j) + vstar_lagF(k) * delta_hxy(xGrid,yGrid,xLag,yLag)*h_l**2*dt
                enddo
            enddo
        enddo

        !call output_force_grid(frame,t,u_star,v_star)

    ! ===================================
        ! Solve projection poisson problem
        call bc(u_star,v_star,U_inf)
        forall(i=1:N_x,j=1:N_y)
            ! page 16 lecture 4 notes
            Q(i,j) = 1.d0/dt * ((u_star(i,j)-u_star(i-1,j)) / h + (v_star(i,j)-v_star(i,j-1)) / h * y_center(j))
        end forall
        ! Solve poisson problem
        call solve_poisson(p,Q,a,b,cm,cp)

    ! ===================================
        ! Update velocities to end time
        ! see work on page 18 lecture 4
        forall (i=1:N_x,j=1:N_y)
            ! appears rho is 1.d0 from Q matrix
            u(i,j) = u_star(i,j) - dt * (p(i+1,j) - p(i,j)) / (rho*h)
            v(i,j) = v_star(i,j) - dt * (p(i,j+1) - p(i,j)) * y_edge(j) / (rho*h)
        end forall

    ! ===================================
        ! Check convergence
        R = 0.d0
        do j=1,N_y
            do i=1,N_x
                if (R < abs(u(i,j)-u_old(i,j)) .or. R < abs(v(i,j)-v_old(i,j))) then
                    R = max(R,abs(u(i,j)-u_old(i,j)),abs(v(i,j)-v_old(i,j)))
                    i_R = i
                    j_R = j
                endif
            enddo
        enddo

        ! Finish up loop
        !print "(a,i4,a,i3,a,i3,a,e16.8)","Loop ",n,": (",i_R,",",j_R,") - R = ",R
        write (13,"(i5,i4,i4,e16.8)") n,i_R,j_R,R
        ! Write out u,v,p every n_steps
        if (mod(n,n_steps) == 0) then
            frame = frame + 1
            call output_grid(frame,t,u,v,p)
            print *, "R = ", R
            print "(a,i3,a,i4,a,e16.8)","Writing frame ",frame," during step n=",n," t=",t
        endif
        ! Check tolerance
        if (R < TOLERANCE) then
            print *, "Convergence reached, R = ",R
            call output_grid(frame,t,u,v,p)
            print "(a,i3,a,i4,a,e16.8)","Writing frame ",frame," during step n=",n," t=",t
            exit
        endif
        ! We did not reach our tolerance, iterate again
        t = t + dt
    enddo
    if (R > TOLERANCE) then
        print "(a,e16.8)","Convergence was never reached, R = ", R

        call output_grid(frame,t,u,v,p) ! ouput last grid?
    endif
    print "(a,i3,a,i4,a,e16.8)","Tolerance reaced!, Writing frame ",frame," during step n=",n," t=",t
    close(13)
    close(70)
    !close(77)
end program main


! ===================================
! Boundary Conditions
! ===================================

subroutine bc(u,v,U_inf)

    use grid_module

    implicit none

    ! Input/Output arguments
    double precision, intent(in) :: U_inf
    double precision, intent(inout) :: u(0:N_x+1,0:N_y+1)
    double precision, intent(inout) :: v(0:N_x+1,0:N_y+1)

    ! Locals
    integer :: i,j

    ! Lower Boundary         Upper Boundary
    !  i,1  i,1    |
    !   x    +     x i,1     |   o   |
    !   |          |         |i,N_y+1|
    ! -------o--------       x   +   x
    !  |    i,0    |         |       |
    !  x     +     x        -|---o---|-
    !  |    i,0   i,0        | i,N_y |
    !  |           |         x   +   x
    !                               i,N_y
    ! Wall                  Free shear
    ! u(i,0) = -u(i,1)      u(i,N_y+1) = u(i,N_y)
    ! v(i,0) = 0.d0         v(i,N_y+1) = v(i,N_y)
    ! p(i,0) = p(i,1)       p(i,N_y+1) = 0.d0
    forall(i=0:N_x+1)
        ! TA's ?
        !u(i,0) = u(i,1)
        !v(i,0) = v(i,1)
        !u(i,N_y+1) = u(i,N_y)
        !v(i,N_y+1) = v(i,N_y)

        ! given
        !u(i,0) = u(i,1)
        !v(i,0) = v(i,1)

        u(i,0) = u(i,1)
        v(i,0) = v(i,1)
        u(i,N_y+1) = u(i,N_y)
        v(i,N_y+1) = v(i,N_y)

    end forall
    ! Left Boundaries
    !   x 0,j  |      x 1,j
    !  0,j     |                        u(0,j) = U_inf
    !   +      o 0,j  + 1,j   o 1,j     v(0,j) = 0.d0
    !          |                        p(0,j) = p(1,j)  (P_x = 0)
    ! Right Boundaries (outflow)
    !     x N_x,j   |         x N_x+1,j             u(N_x+1,j) = 2 u(N_x,j) - u(N_x-1,j)
    !               |                               v(N_x+1,j) = 2 v(N_x,j) - v(N_x-1,j)
    !     + N_x,j   o N_x,j   + N_x+1,j  o N_x+1,j  p(N_x+1,j) = p(N_x,j)
    !               |
    forall(j=0:N_y+1)
        ! TAs
        !u(0,j) = U_inf
        !v(0,j) = 0.d0
        !u(N_x+1,j) = u(N_x,j)
        !v(N_x+1,j) = v(N_x,j)

        ! given
        u(0,j) = U_inf
        v(0,j) = 0.d0

        u(N_x+1,j) = 2.d0*u(N_x,j) - u(N_x-1,j)
        !u(N_x+1,j) = u(N_x,j)
        !v(N_x+1,j) = v(N_x,j)
        v(N_x+1,j) = 2.d0*v(N_x,j) - v(N_x-1,j)
    end forall

end subroutine bc


! ===================================
! Solve the poisson problem
!  laplace(P)_ij = Q_ij
! ===================================
subroutine solve_poisson(P,Q,a,b,cm,cp)

    use grid_module

    implicit none

    ! Input
    double precision, intent(in) :: Q(1:N_x,1:N_y),a
    double precision, intent(in) :: b(1:N_y),cm(1:N_y),cp(1:N_y)
    double precision, intent(inout) :: P(0:N_x+1,0:N_y+1)

    ! Solver parameters
    logical, parameter :: verbose = .false.
    integer, parameter :: MAX_ITERATIONS = 10000
    double precision, parameter :: TOLERANCE = 10.d-4
    double precision, parameter :: w = 1.6d0 ! 1 (GS) < w < 2

    ! Local variables
    integer :: i,j,n
    double precision :: R,P_old

    do n=1,MAX_ITERATIONS
        R = 0.d0
        ! Boundary conditions
        forall (j=0:N_y+1) ! left and righ
            ! TAs
            !P(0,j) = P(1,j)
            !P(N_x+1,j) = 0

            ! given
            P(0,j) = P(1,j)        ! Left
            !P(N_x+1,j) = P(N_x,j)  ! Right

            ! iman!
            P(N_x+1,j) = P(N_x, j)
        end forall
        forall (i=0:N_x+1) ! top and bottom
            !TAs
            !P(i,0) = P(i,1)
            !P(i,N_y+1) = P(i,N_y)

            ! given
            !P(i,0) = P(i,1)        ! Bottom wall
            P(i,0) = 0.d0          ! Free stream
            P(i,N_y+1) = 0.d0      ! Free stream

        end forall

        do j=1,N_y
            do i=1,N_x
                ! Update p
                P_old = P(i,j)
                P(i,j) = (Q(i,j) - (a * (P(i+1,j) + P(i-1,j)) + cp(j)*P(i,j+1) + cm(j)*P(i,j-1))) / b(j)
                ! relaxation
                P(i,j) = w * P(i,j) + (1-w) * P_old
                ! Calculate Residual
                R = max(R,abs(P(i,j)-P_old))
            enddo
        enddo
        !print *, "R: ", R
        ! Print out convergence and exit
        if (verbose) then
            print "(a,i5,a,e16.8,a)", "(",n,",",R,")"
        endif
        if (R < TOLERANCE) then
            exit
        endif
    enddo
    ! Check to see if we converged and quit if we did not
    if (n > MAX_ITERATIONS) then
        print *,"Poisson solver did not converge, exiting!"
        print "(a,e16.8)","R = ",R
        print *,"iteration: ", n
        stop
    else
        !print "(a,i6,a,e16.8)", "Solver converged in ",n," steps: R = ",R
        ! Boundary conditions
        ! Boundary conditions
        forall (j=0:N_y+1) ! left and righ
            ! TAs
            !P(0,j) = P(1,j)
            !P(N_x+1,j) = 0

            ! given
            P(0,j) = P(1,j)        ! Left
            !P(N_x+1,j) = P(N_x,j)  ! Right

            ! iman!
            P(N_x+1,j) = P(N_x, j)
        end forall
        forall (i=0:N_x+1) ! top and bottom
            !TAs
            !P(i,0) = P(i,1)
            !P(i,N_y+1) = P(i,N_y)

            ! given
            !P(i,0) = P(i,1)        ! Bottom wall
            P(i,0) = 0.d0          ! Free stream
            P(i,N_y+1) = 0.d0      ! Free stream

        end forall


    !    forall (j=0:N_y+1)
    !        P(0,j) = P(1,j)        ! Left
    !        P(N_x+1,j) = P(N_x,j)  ! Right
    !    end forall
    !    forall (i=0:N_x+1)
    !        P(i,0) = 0.d0        ! Free stream
    !        P(i,N_y+1) = 0.d0      ! Free stream
    !    end forall
    endif

end subroutine solve_poisson

