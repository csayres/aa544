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
    integer :: N_y , N_x, N_leadIn, nlx, nly, n_lagrangian_points
    double precision :: L_x, h, L_y, HL_y, HL_x, H_xOffset, h_l, Re
    character*20 :: fileSuffix

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
    double precision, allocatable :: x_edge(:),x_center(:), x_lag(:), x_euler_per_lag(:,:)
    double precision, allocatable :: y_edge(:),y_center(:), y_lag(:), y_euler_per_lag(:,:)

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

        x_front = H_xOffset - HL_x
        x_back = H_xOffset
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

        ! next for each lagrangian_points find the associated eulearian grid points (x and y separately)

    end subroutine setup_lagrangian_points

    ! Setup and allocate all lagrangian points
    subroutine setup_lagrangian_points_circle()

        implicit none


        !locals
        double precision :: x_cent, y_cent, r
        double precision :: y_top, y_bottom
        integer :: i, running_i ! number of points in the x and y directions

        !h_l = h/2 ! spacing between lagrangian points

        ! determine number of points in x and y directions
        !nlx = HL_x / h_l
        !nly = HL_y / h_l

        ! all lagrangian points
        allocate(x_lag(1:(n_lagrangian_points)))
        allocate(y_lag(1:(n_lagrangian_points)))

        x_cent = H_xOffset + 0.5*HL_x ! x center
        y_cent = 0.5 * L_y ! y center
        r = 0.5*HL_y ! radius

        do i=1,n_lagrangian_points
            x_lag(i) = r*cos(i*1.d0/(n_lagrangian_points*1.d0)*2*3.14159) + x_cent
            y_lag(i) = r*sin(i*1.d0/(n_lagrangian_points*1.d0)*2*3.14159) + y_cent
        enddo


        ! next for each lagrangian_points find the associated eulearian grid points (x and y separately)

    end subroutine setup_lagrangian_points_circle

    subroutine output_lagrangian_points()

        ! local
        integer :: i
        integer :: x_ind, y_ind
        open(unit=35,file='_output/lagrangian_points'//fileSuffix//'.dat',access='sequential',status='unknown')
        do i=1,nlx*2+nly*2
            write(35,*) x_lag(i), y_lag(i)
        enddo
        close(35)

    end subroutine output_lagrangian_points

    subroutine setup_subDomain()

        ! local
        integer :: i
        integer :: x_ind, y_ind

        ! four indices per lagrangian point (in x and y)
        allocate(x_euler_per_lag(1:n_lagrangian_points, 1:5))
        allocate(y_euler_per_lag(1:n_lagrangian_points, 1:5))
        do i=1,n_lagrangian_points
            ! get the closest indices to lagrangian point on the eulerian grid
            x_ind = x_lag(i)/h
            y_ind = y_lag(i)/h
            ! add 2 on either side
            y_euler_per_lag(i,1) = y_ind - 2
            y_euler_per_lag(i,2) = y_ind - 1
            y_euler_per_lag(i,3) = y_ind
            y_euler_per_lag(i,4) = y_ind + 1
            y_euler_per_lag(i,5) = y_ind + 2

            x_euler_per_lag(i,1) = x_ind - 2
            x_euler_per_lag(i,2) = x_ind - 1
            x_euler_per_lag(i,3) = x_ind
            x_euler_per_lag(i,4) = x_ind + 1
            x_euler_per_lag(i,5) = x_ind + 2
        end do

    end subroutine setup_subDomain

    subroutine output_grid_centers()

        ! local
        integer :: i
        open(unit=35,file='_output/x_points'//fileSuffix//'.dat',access='sequential',status='unknown')
        do i=1,N_x
            write(35,*) x_center(i)
        enddo
        close(35)
        open(unit=35,file='_output/y_points'//fileSuffix//'.dat',access='sequential',status='unknown')
        do i=1,N_y
            write(35,*) y_center(i)
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
        write(70,*) "Grid Size", N_x, N_y, "Reynolds = ", Re, "height= ", HL_y, "width= ", HL_x
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
        write(77,*) "Grid Size", N_y, N_x, "Reynolds = ", Re, "height= ", HL_y, "width= ", HL_x
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

!calculate values for a tensor
subroutine geta_11(i,j,nx,ny,h,u,v,nu_t,a_11)
    implicit none
    integer, intent(in) :: i,j,nx,ny
    double precision, intent(in) :: h
    double precision, dimension(0:nx+1,0:ny+1), intent(in) :: u,v,nu_t
    double precision, intent(out)::a_11
    a_11 = -2.d0*nu_t(i,j)*(u(i,j)-u(i-1,j))/h
end subroutine geta_11

subroutine geta_22(i,j,nx,ny,h,u,v,nu_t,a_22)
    implicit none
    integer, intent(in) :: i,j,nx,ny
    double precision, intent(in) :: h
    double precision, dimension(0:nx+1,0:ny+1), intent(in) :: u,v,nu_t
    double precision, intent(out)::a_22
    a_22 = -2.d0*nu_t(i,j)*(v(i,j)-v(i,j-1))/h
end subroutine geta_22

subroutine geta_12(i,j,nx,ny,h,u,v,nu_t,a_12)
    implicit none
    integer, intent(in) :: i,j,nx,ny
    double precision, intent(in) :: h
    double precision, dimension(0:nx+1,0:ny+1), intent(in) :: u,v,nu_t
    double precision, intent(out)::a_12
    a_12 = -(nu_t(i,j) + nu_t(i+1,j) + nu_t(i,j+1) + nu_t(i+1,j+1))/4.d0*(u(i,j+1)-u(i,j) + v(i+1,j)-v(i,j))/h
end subroutine geta_12

subroutine geta_21(i,j,nx,ny,h,u,v,nu_t,a_21)
    implicit none
    integer, intent(in) :: i,j,nx,ny
    double precision, intent(in) :: h
    double precision, dimension(0:nx+1,0:ny+1), intent(in) :: u,v,nu_t
    double precision, intent(out)::a_21
    call geta_12(i,j,nx,ny,h,u,v,nu_t,a_21)
end subroutine geta_21

subroutine nu_t_in(i,j,nx,ny, h, y,P,u,v,rho,nu,nu_t)
    implicit none
    integer, intent(in) :: i,j,nx,ny
    double precision, intent(in) :: h, rho, nu
    double precision, dimension(0:ny+1), intent(in) :: y
    double precision, dimension(0:nx+1,0:ny+1), intent(in) :: u,v,P
    double precision, intent(out)::nu_t
    !local
    double precision :: du_dy, dv_dx, l
    du_dy = ((u(i,j+1) + u(i-1,j-1))/2.d0 - (u(i,j-1) + u(i-1,j-1))/2.d0 )/(2.d0*h)
    dv_dx = ((v(i-1,j-1)+v(i-1,j))/2.d0 - (v(i+1,j)+v(i+1,j-1))/2.d0 )/(2.d0*h)
    call lm(i,j,y,P,u,rho, nu, l)
    nu_t = l**2*sqrt(du_dy**2+dv_dx**2)
end subroutine nu_t_in

subroutine getTao_omega(i,nx,ny,h,u,rho, nu,tao_omega) !wall sheer stress. evaluated at y=0
    implicit none
    integer, intent(in) :: i,nx,ny
    double precision, intent(in) :: h,rho, nu
    double precision, dimension(0:nx+1,0:ny+1)::u
    double precision, intent(out) :: tao_omega
    !local
    double precision :: du_dy
    du_dy = ((u(i,1) + u(i-1,1))/2.d0 - (u(i,0) + u(i-1,0))/2.d0 )/(h)

    tao_omega = rho*nu*du_dy
end subroutine getTao_omega

subroutine getU_tao(i,nx,ny,h,u,rho, nu,u_tao) !friction velocity
    implicit none
    integer, intent(in) :: i,nx,ny
    double precision, intent(in) :: h,rho, nu
    double precision, dimension(0:nx+1,0:ny+1)::u
    double precision, intent(out)::u_tao
    !local
    double precision :: t_w
    call getTao_omega(i,nx,ny,h,u,rho, nu,t_w)
    u_tao = sqrt(t_w/rho)
end subroutine getU_tao

subroutine getA_plus(i,j,nx,ny,h,y,P,u,rho,nu,a_plus) !friction velocity
    implicit none
    integer, intent(in) :: i,j,nx,ny
    double precision, intent(in) :: h, rho, nu
    double precision, dimension(0:ny+1), intent(in) :: y
    double precision, dimension(0:nx+1,0:ny+1), intent(in) :: u,P
    double precision, intent(out) :: a_plus
    !local
    double precision :: dp_dx, u_t,t_w
    call getU_tao(i,nx,ny,h,u,rho, nu, u_t)
    dp_dx = (P(i+1,j)-P(i-1,j))/(2.0*h)
    call getTao_omega(i,nx,ny,h,u,rho, nu,t_w)
    a_plus = 26.d0/sqrt(1.d0+y(j)*dp_dx/(rho*t_w**2))
end subroutine getA_plus

subroutine getY_plus(i,j,nx,ny,h,y,u,rho, nu,y_plus) !friction velocity
    implicit none
    integer, intent(in) :: i,j,nx,ny
    double precision, intent(in) :: h, rho, nu
    double precision, dimension(0:ny+1), intent(in) :: y
    double precision, dimension(0:nx+1,0:ny+1), intent(in) :: u
    double precision, intent(out) :: y_plus
    !local
    double precision:: t_w
    call getTao_omega(i,nx,ny,h,u,rho, nu,t_w)
    y_plus = y(j) / (nu*sqrt(rho/t_w))
end subroutine getY_plus

subroutine lm(i,j,nx,ny,h,y,P,u,rho, nu, l) !friction velocity
    implicit none
    integer, intent(in) :: i,j,nx,ny
    double precision, intent(in) :: h, rho, nu
    double precision, dimension(0:ny+1), intent(in) :: y
    double precision, dimension(0:nx+1,0:ny+1), intent(in) :: u,P
    double precision, intent(out) :: l
    !local
    double precision :: K, yp, ap
    K = 0.40
    call getY_plus(i,j,nx,ny,h,y,u,rho, nu,yp)
    call getA_plus(i,j,nx,ny,h,y,P,u,rho,nu, ap)
    l = K*y(j)*(1.d0 - exp(-yp/ap))
end subroutine lm

subroutine nu_t_out(i,j,nx,ny,h,y,P,u,v,rho,nu,Ue, nu_t)
    implicit none
    integer, intent(in) :: i,j,nx,ny
    double precision, intent(in) :: h, rho, nu, Ue
    double precision, dimension(0:ny+1), intent(in) :: y
    double precision, dimension(0:nx+1,0:ny+1), intent(in) :: u,P,v
    double precision, intent(out) :: nu_t
    !local
    double precision :: alpha, mt, fk
    alpha = 0.0168
    call momentumThick(i,nx,ny,h,u,Ue, mt)
    call fkleb(i,j,nx,ny,y,u,Ue,fk)
    nu_t = alpha*Ue*mt*fk
end subroutine nu_t_out

subroutine momentumThick(i,nx,ny,h,u,Ue, mt)
    implicit none
    integer, intent(in) :: i,nx,ny
    double precision, intent(in) :: h, Ue
    double precision, dimension(0:nx+1,0:ny+1), intent(in) :: u
    double precision, intent(out) :: mt
    !local
    double precision :: theSum
    integer :: j
    theSum = 0.d0
    ! intetrate over all at center location i (average left and right u)
    do j=1,ny
        theSum = theSum + (1 - (u(i,j)+u(i-1,j))/2.0/Ue)*h
    enddo
    mt = theSum
end subroutine momentumThick

subroutine flkeb(i,j,nx,ny,y,u,Ue,fk)
    implicit none
    integer, intent(in) :: i,j,nx,ny
    double precision, intent(in) :: Ue
    double precision, dimension(0:ny+1), intent(in) :: y
    double precision, dimension(0:nx+1,0:ny+1), intent(in) :: u
    double precision, intent(out) :: fk
    !local
    double precision :: delta
    integer :: jj
    ! intetrate over all at center location i (average left and right u), find out where it approximately
    ! equals Ue
    do jj=1,ny
        delta = (u(i,jj)+u(i-1,jj))/2.0
        if (delta >= 0.99*Ue) exit
    enddo
    ! delta is a y value
    delta = y(jj)
    fk = 1.d0/(1.d0 + 5.5 * (y(j)/delta)**6)
end subroutine flkeb


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
    integer, parameter :: MAX_ITERATIONS = 10000
    double precision, parameter :: TOLERANCE = 1d-4, CFL = 0.02
    logical, parameter :: write_star = .false.
    integer :: n_steps

    ! ===================================
    ! Physics
    double precision :: U_inf = 1.d0
    double precision :: rho,nu !rho = 1.d0, nu=1.d-3

    ! ===================================
    ! Velocity and pressures
    double precision, allocatable :: u(:,:),v(:,:),p(:,:),u_star(:,:),v_star(:,:)
    double precision, allocatable :: u_old(:,:), v_old(:,:)
    integer, allocatable :: nu_t_flag(:) ! boolean array, for using nu_t in or out 0 = in 1 = out

    ! ===================================
    ! Locals
    character*20 :: arg
    integer :: i,ii,j,jj,n,m,frame,i_R,j_R, k, boxLen, maxIterNorm
    double precision :: R,t,dt,a, lagFuSum, lagFvSum
    double precision, allocatable :: Flux_ux(:,:),Flux_uy(:,:),Flux_vy(:,:),nu_t(:,:)
    double precision :: uu_x,uv_y,uv_x,vv_y,u_xx,u_yy,v_xx,v_yy
    double precision, allocatable :: Q(:,:),b(:),cp(:),cm(:), ustar_lagF(:), vstar_lagF(:)
    double precision :: xGrid, yGrid, xLag, yLag, turbViscIn, turbViscOut, turbVisc
    double precision :: a_11a,a_12a,a_21a,a_22a,a_11b,a_12b,a_21b,a_22b
    double precision :: da11_dx, da12_dy, da21_dx, da22_dy
    ! ===================================

    ! Get command line arguments
    if (iargc() /= 3) then
        print *,"Wrong number of command line arguments, expected 2."
        print *,"    Lx - Length of domain in x"
        print *,"    boxLen - Length of box"
        print *," Re = reynolds number"
        print *, "got ", iargc()
        stop
    else
        call getarg(1,arg)
        read (arg,'(I10)') N_y
        call getarg(2,arg)
        read (arg,'(I10)') boxLen
        call getarg(3,arg)
        read (arg,'(F7.0)') Re
    endif
    write(fileSuffix, "(I3,'_',I1,'_',F7.0)") N_y, boxLen, Re
    print *, fileSuffix
    print *, Re

!!!!!!!!!!!!!!!!!!
!!  Parameters: !!
!!!!!!!!!!!!!!!!!!

    !N_x=10  !Number of grid points in x-direction
    !N_y = 128   !Number of grid points in y-direction
    L_x = 300 !Length of box in x-direction
    L_y = 300  !Length of box in y-direction


    n_steps = MAX_ITERATIONS/100 !Interval that u,v and p are printed to UVP.dat




    ! Setup grid and arrays
    call setup_grid()
    call output_grid_centers()

    !!! Lagrangian Points
    HL_y = 10 !0.05 * L_y  ! Length of rect bluff y-direction
    HL_x = boxLen*HL_y ! Length of rect bluff in x-direction
    H_xOffset = 0.5 * L_x ! how far along x before bluff starts, bluff will always be
                           ! centered in the domain.
    h_l = h ! spacing of lagrangian points
    ! determine number of points in x and y directions
    nlx = HL_x / h_l
    nly = HL_y / h_l
    n_lagrangian_points = 2*nlx + 2*nly
    !call setup_lagrangian_points()
    call setup_lagrangian_points_circle()
    call output_lagrangian_points()
    call setup_subDomain()
    N_leadIn = N_x / 25
    ! ===================================

    allocate(Flux_ux(1:N_x+1,0:N_y+1))
    allocate(Flux_uy(0:N_x,0:N_y))
    allocate(Flux_vy(0:N_x+1,0:N_y))
    allocate(Q(1:N_x,1:N_y))
    allocate(b(1:N_y),cp(1:N_y),cm(1:N_y))

    ! Calculate matrix coefficients for Poisson solve
    a = 1.d0 / h**2
    forall (j=1:N_y)
        b(j) = - (2.d0 / h**2 + y_center(j) / h**2 * (1.d0 + 1.d0))
        cp(j) = (1.d0 * 1.d0) / h**2
        cm(j) = (1.d0 * 1.d0) / h**2
    end forall

    ! Velocity and pressure arrays
    allocate(u(0:N_x+1,0:N_y+1),u_star(0:N_x+1,0:N_y+1))
    allocate(v(0:N_x+1,0:N_y+1),v_star(0:N_x+1,0:N_y+1))
    allocate(p(0:N_x+1,0:N_y+1))
    allocate(nu_t(0:N_x+1,0:N_y+1))
    allocate(u_old(1:N_x,1:N_y),v_old(1:N_x,1:N_y))
    allocate(nu_t_flag(1:N_x))
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
    nu_t = 0.d0

    !Re = 100.0
    nu = (U_inf*HL_y)/Re
    !nu = 1.d0/Re
    rho = 1.d0
    !dt =  CFL * h / 500.d0
    !dt = CFL * h / (Re * U_inf)
    dt = CFL * h / (U_inf)
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
        nu_t_flag = 0
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
                ! turbulent viscocity terms
                ! first determine whether to use nu_t in or out
                call nu_t_out(i,j,N_x,N_y,h,y_center,P,u,v,rho,nu,U_inf,turbVisc)
                if (nu_t_flag(j)==0) then
                    ! calculate both nu_t versions
                    call nu_t_in(i,j,N_x,N_y, h, y_center,P,u,v,rho,nu,turbViscIn)
                    if (abs(1.d0 - turbVisc/turbViscIn)<=0.1) then
                        ! switch to outer
                        nu_t_flag(j)=1
                    else
                        turbVisc = turbViscIn
                    end if
                end if
                nu_t(i,j) = turbVisc
            enddo
        enddo
        do j=1,N_y
            do i=1,N_x
                ! Advective terms, see Lecture 4 notes, page 12
                uu_x = (Flux_ux(i+1,j) - Flux_ux(i,j)) / h
                uv_y = 1.d0 * (Flux_uy(i,j) - Flux_uy(i,j-1)) / h !

                uv_x = (Flux_uy(i,j) - Flux_uy(i-1,j)) / h
                vv_y = 1.d0 * (Flux_vy(i,j) - Flux_vy(i,j-1)) / h

                ! Diffusive terms
                u_xx = (u(i+1,j) - 2.d0*u(i,j) + u(i-1,j)) / (h**2)
                u_yy = (1.d0 / h**2) * (1.d0 * (u(i,j+1) - u(i,j)) - 1.d0 * (u(i,j) - u(i,j-1))) !edge leads center by 0.5, verify by looking at poisson solver
                v_xx = (v(i+1,j) - 2.d0*v(i,j) + v(i-1,j)) / (h**2)
                v_yy = (1.d0 / h**2) * (1.d0 * (v(i,j+1) - v(i,j)) - 1.d0 * (v(i,j) - v(i,j-1)))

                !! solve for a_ij's
                call geta_11(i+1,j,N_x,N_y,h,u,v,nu_t,a_11a)
                call geta_11(i-1,j,N_x,N_y,h,u,v,nu_t,a_11b)

                call geta_12(i,j+1,N_x,N_y,h,u,v,nu_t,a_12a)
                call geta_12(i,j-1,N_x,N_y,h,u,v,nu_t,a_12b)

                call geta_21(i+1,j,N_x,N_y,h,u,v,nu_t,a_21a)
                call geta_21(i-1,j,N_x,N_y,h,u,v,nu_t,a_21b)

                call geta_22(i,j+1,N_x,N_y,h,u,v,nu_t,a_22a)
                call geta_22(i,j-1,N_x,N_y,h,u,v,nu_t,a_22b)

                ! solve d_a11/dx
                da11_dx = (a_11a - a_11b)/(2*h)
                da12_dy = (a_12a - a_12b)/(2*h)
                da21_dx = (a_21a - a_21b)/(2*h)
                da22_dy = (a_22a - a_22b)/(2*h)

                ! Update to u* and v* values
                u_star(i,j) = u(i,j) + dt * (-(uu_x + uv_y + da11_dx + da12_dy) + nu*(u_xx + u_yy))
                v_star(i,j) = v(i,j) + dt * (-(uv_x+vv_y + da21_dx + da22_dy) + nu*(v_xx + v_yy))
            enddo
        enddo

        ! Debug, save out u_star,v_star,p
        if (write_star) then
            frame = frame + 1
            print "(a,i3,a,i4,a,e16.8)","Writing frame ",frame," during step n=",n," t=",t
            call output_grid(frame,t,u_star,v_star,p)
        endif
    !===================================

    ! ===================================
        ! Solve projection poisson problem
        call bc(u_star,v_star,U_inf)
        forall(i=1:N_x,j=1:N_y)
            ! page 16 lecture 4 notes
            Q(i,j) = 1.d0/dt * ((u_star(i,j)-u_star(i-1,j)) / h + (v_star(i,j)-v_star(i,j-1)) / h * 1.d0)
        end forall
        ! Solve poisson problem
        call solve_poisson(p,Q,a,b,cm,cp)

    ! ===================================
        ! Update velocities to end time
        ! see work on page 18 lecture 4
        forall (i=1:N_x,j=1:N_y)
            ! appears rho is 1.d0 from Q matrix
            u(i,j) = u_star(i,j) - dt * (p(i+1,j) - p(i,j)) / (rho*h)
            v(i,j) = v_star(i,j) - dt * (p(i,j+1) - p(i,j)) * 1.d0 / (rho*h)
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

        ! Write out u,v,p every n_steps
        if (mod(n,n_steps) == 0) then
            frame = frame + 1
            call output_grid(frame,t,u,v,p)
            print *, "R = ", R
            write (13,"(i5,i4,i4,e16.8)") n,i_R,j_R,R
            print "(a,i3,a,i4,a,e16.8)","Writing frame ",frame," during step n=",n," t=",t
        endif
        ! Check tolerance
        !if (R < TOLERANCE) then
            !print *, "Convergence reached, R = ",R
            !call output_grid(frame,t,u,v,p)
            !print "(a,i3,a,i4,a,e16.8)","Writing frame ",frame," during step n=",n," t=",t
            !exit
        !endif
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
    call output_grid(frame,t,u,v,p)
    print *, "!!!!!!!!!!!!!!!!!!!1END!!!!!!!!!!!!!!!", fileSuffix
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !u(i,0) = u(i,1)
        !v(i,0) = v(i,1)
        u(i,0) = -u(i,1)
        v(i,0) = 0.d0
        u(i,N_y+1) = u(i,N_y)
        v(i,N_y+1) = v(i,N_y)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end forall
 !   forall(i=0:N_leadIn)
 !       u(i,0) = u(i,1)
 !       v(i,0) = v(i,1)
 !   end forall
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
    double precision, parameter :: TOLERANCE = 100.d-4
    double precision, parameter :: w = 1.6d0 ! 1 (GS) < w < 2

    ! Local variables
    integer :: i,j,n
    double precision :: R,P_old

    do n=1,MAX_ITERATIONS
        R = 0.d0
        ! Boundary conditions

        forall (i=0:N_x+1) ! top and bottom
            !TAs
            !P(i,0) = P(i,1)
            !P(i,N_y+1) = P(i,N_y)

            ! given
            P(i,0) = P(i,1)        ! Bottom wall
            !P(i,0) = 0.d0          ! Free stream
            P(i,N_y+1) = 0.d0      ! Free stream


        end forall
    !    forall(i=0:N_leadIn)
    !        P(i,0) = 0.d0
    !    end forall
        forall (j=0:N_y+1) ! left and righ
            ! TAs
            !P(0,j) = P(1,j)
            P(N_x+1,j) = 0

            ! given
            P(0,j) = P(1,j)        ! Left
            !P(N_x+1,j) = P(N_x,j)  ! Right

            ! iman!
            !P(N_x+1,j) = P(N_x, j)
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
        forall (i=0:N_x+1) ! top and bottom
            !TAs
            !P(i,0) = P(i,1)
            !P(i,N_y+1) = P(i,N_y)

            ! given
            P(i,0) = P(i,1)        ! Bottom wall
            !P(i,0) = 0.d0          ! Free stream
            P(i,N_y+1) = 0.d0      ! Free stream


        end forall
    !    forall(i=0:N_leadIn)
    !        P(i,0) = 0.d0
    !    end forall
        forall (j=0:N_y+1) ! left and righ
            ! TAs
            !P(0,j) = P(1,j)
            P(N_x+1,j) = 0

            ! given
            P(0,j) = P(1,j)        ! Left
            !P(N_x+1,j) = P(N_x,j)  ! Right

            ! iman!
            !P(N_x+1,j) = P(N_x, j)
        end forall
    endif

end subroutine solve_poisson

