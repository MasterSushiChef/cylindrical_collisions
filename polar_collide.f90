! ------------------ FUNCTIONS AND SUBROUTINES ------------------ !
module polar_collide

external :: dgesv
double precision, parameter :: Pi = 3.14159265358979323d0

contains

! ------------------ COLLIDE SUBROUTINES ------------------ !
subroutine build_cdf(n_r, n_theta, vdf, cdf)
    implicit none
    double precision, allocatable, intent(in) :: vdf(:,:)
    integer, intent(in) :: n_r, n_theta
    double precision, allocatable, intent(inout) :: cdf(:)

    integer :: i, j, index

    ! Build cumulative distribution function.
    cdf(1) = vdf(1, 1)
    index = 2

    do j = 1,n_theta
        do i = 2,n_r
            cdf(index) = cdf(index-1) + vdf(i, j)
            index = index + 1
        end do
    end do
end subroutine build_cdf

! Deplete mass from vdf and calculate pre-collision velocity index.
subroutine precollision(grid_r, grid_theta, n_r, n_theta, cdf, vr, vtheta)
    implicit none
    double precision, allocatable, intent(in) :: grid_r(:), grid_theta(:)
    integer, intent(in) :: n_r, n_theta ! size of grid
    double precision, allocatable, intent(in) :: cdf(:)
    double precision, intent(out) :: vr, vtheta

    double precision :: random_num, max_cdf, find ! variables for finding random index
    integer :: l ! location of velocity
    integer :: x, lower, upper
    integer :: i, j ! index to r, theta

    ! Calculate value to search for in cdf to determine vr, thetav.
    call random_init(.false., .true.)
    call random_number(random_num)
    max_cdf = cdf(ubound(cdf,1))
    find = max_cdf * random_num

    ! Find velocity values from index in cdf and deplete delta_m.
    l = binary_search(cdf, find)

    if (l .eq. 1) then
        i = 1
        j = 1
    else
        do x = 0, n_theta-1
            lower = 2 + x * (n_r - 1)
            upper = 2 + (x + 1) * (n_r - 1)
            if (l .ge. lower .and. l .lt. upper) then
                j = x+1
                i = l - lower + 2
            end if
        end do
    end if

    vr = grid_r(i)
    vtheta = grid_theta(j)
end subroutine precollision

! Given two points, calculate and return post collision velocity.
subroutine collide(vr1, vr2, vtheta1, vtheta2, vr1_prime, vr2_prime, vtheta1_prime, vtheta2_prime)
    implicit none
    double precision, intent(in) :: vr1, vr2, vtheta1, vtheta2 ! pre collision velocity index
    double precision, intent(out) :: vr1_prime, vr2_prime, vtheta1_prime, vtheta2_prime ! post colision velocity index

    double precision :: vx1, vy1, vx2, vy2 ! pre collision velocities
    double precision :: vx1_prime, vy1_prime, vx2_prime, vy2_prime ! post collision velocities
    double precision :: g ! relative velocity magnitude
    double precision :: wx, wy ! center of mass velocity
    double precision :: gx_prime, gy_prime
    double precision :: rf, phi, cos_theta, sin_theta ! collision sphere rotation angles
    double precision :: delta_m = 1.3503331271824742D-006 ! arbitrary
    double precision :: p_x1, p_y1, p_x2, p_y2, e_1, e_2
    
    ! Precalculate x, y velocities.
    vx1 = vr1 * cos(vtheta1)
    vy1 = vr1 * sin(vtheta1)
    vx2 = vr2 * cos(vtheta2)
    vy2 = vr2 * sin(vtheta2)
    if (abs(vx1) .lt. 1D-6) vx1 = 0.0
    if (abs(vy1) .lt. 1D-6) vy1 = 0.0
    if (abs(vx2) .lt. 1D-6) vx2 = 0.0
    if (abs(vy2) .lt. 1D-6) vy2 = 0.0

    p_x1 = delta_m * (vx1 + vx2)
    p_y1 = delta_m * (vy1 + vy2)
    e_1 = vx1**2 + vy1**2 + vx2**2 + vy2**2

    ! Calculate relative velocity magnitude and center of mass velocity of collision pair.
    g = sqrt((vx1 - vx2)**2 + (vy1 - vy2)**2)
    wx = 0.5*(vx1 + vx2)
    wy = 0.5*(vy1 + vy2)

    ! Calculate scattering angle.
    call random_init(.false., .true.)
    call random_number(rf)
    phi = 2*Pi*rf ! TODO: some angles not necessary
    cos_theta = 2*rf - 1
    sin_theta = sqrt(1 - cos_theta**2)

    ! Calculate post coliision velocities in x,y.
    gx_prime = 0.5 * g*sin_theta*cos(phi) ! TODO: change these two things
    gy_prime = 0.5 * g*sin_theta*sin(phi)

    vx1_prime = wx + gx_prime
    vx2_prime = wx - gx_prime
    vy1_prime = wy + gy_prime
    vy2_prime = wy - gy_prime

    p_x2 = delta_m * (vx1_prime + vx2_prime)
    p_y2 = delta_m * (vy1_prime + vy2_prime)
    e_2 = vx1_prime**2 + vy1_prime**2 + vx2_prime**2 + vy2_prime**2

    if (abs(p_x2 - p_x1) .gt. 1D-12) print *, p_x1, p_x2
    if (abs(p_y2 - p_y1) .gt. 1D-12) print *, p_y1, p_y2
    if (abs(e_2 - e_1) .gt. 1D-12) print *, "energy not conserved"

    ! Convert Cartesian to polar coordinates.
    vtheta1_prime = atan2(vy1_prime, vx1_prime)
    vtheta2_prime = atan2(vy2_prime, vx2_prime)
    vr1_prime = sqrt(vx1_prime**2 + vy1_prime**2)
    vr2_prime = sqrt(vx2_prime**2 + vy2_prime**2)

    ! Turn negative angles into positive ones.
    if (vtheta1_prime .lt. 0) vtheta1_prime = vtheta1_prime + 2*Pi
    if (vtheta2_prime .lt. 0) vtheta2_prime = vtheta2_prime + 2*Pi
end subroutine collide

! Given post collision velocity, calculate which points to map mass back to.
subroutine find_points(vr_prime, vtheta_prime, grid_r, grid_theta, map_coords)
    implicit none
    double precision, allocatable, intent(in) :: grid_r(:), grid_theta(:)
    double precision, intent(in) :: vr_prime, vtheta_prime
    double precision, intent(out) :: map_coords(4,2) ! points to map back to

    double precision :: vtheta_bounds(2), vr_bounds(2) ! grid values bounding v_prime components
    integer :: theta_loc, i

    vtheta_bounds = find_theta_bounds(vtheta_prime, grid_theta)

    ! Determine closest grid points to remap to.
    if (vr_prime .lt. 1.40*grid_r(2)/2) then
        ! If the point is closer to lower vz bound.
        map_coords(1,:) = (/ grid_r(2), vtheta_bounds(1) /)
        map_coords(2,:) = (/ grid_r(2), vtheta_bounds(2) /)

        if (vtheta_bounds(2) .eq. 0.0d0) vtheta_bounds(2) = 2*Pi ! Set bound to 2*pi temporarily to ensure check works.
        if (vtheta_prime - vtheta_bounds(1) .lt. vtheta_bounds(2) - vtheta_prime) then
            ! If the point is closer to lower grid_theta bound.
            if (abs(vtheta_bounds(2) - 2*Pi) .lt. 1D-8) vtheta_bounds(2) = 0.0d0

            theta_loc = find_loc(grid_theta, vtheta_bounds(1)) - 20
            if (theta_loc .lt. 1) then
                theta_loc = size(grid_theta) + theta_loc
            end if
            
            map_coords(3,:) = (/ grid_r(1), 0.0d0 /)
            map_coords(4,:) = (/ grid_r(2), grid_theta(theta_loc) /)
        else
            ! If the point is closer to higher grid_theta bound.
            if (abs(vtheta_bounds(2) - 2*Pi) .lt. 1D-8) vtheta_bounds(2) = 0.0d0

            theta_loc = find_loc(grid_theta, vtheta_bounds(2)) + 20
            if (theta_loc .gt. size(grid_theta)) then
                theta_loc = theta_loc - size(grid_theta)
            end if
            
            map_coords(3,:) = (/ grid_r(1), 0.0d0 /)
            map_coords(4,:) = (/ grid_r(2), grid_theta(theta_loc) /)
        end if
    else if (vr_prime .gt. 1.55*grid_r(2)/2 .and. vr_prime .lt. (grid_r(2) + grid_r(3))/2.0d0) then
        ! If the point is closer to lower vz bound.
        map_coords(1,:) = (/ grid_r(2), vtheta_bounds(1) /)
        map_coords(2,:) = (/ grid_r(2), vtheta_bounds(2) /)
        map_coords(3,:) = (/ grid_r(1), 0.0d0 /)

        if (vtheta_bounds(2) .eq. 0.0d0) vtheta_bounds(2) = 2*Pi ! Set bound to 2*pi temporarily to ensure check works.
        if (vtheta_prime - vtheta_bounds(1) .lt. vtheta_bounds(2) - vtheta_prime) then
            ! If the point is closer to lower grid_theta bound.
            if (abs(vtheta_bounds(2) - 2*Pi) .lt. 1D-8) vtheta_bounds(2) = 0.0d0

            map_coords(4,:) = (/ grid_r(3), vtheta_bounds(1) /)
        else
            ! If the point is closer to higher grid_theta bound.
            if (abs(vtheta_bounds(2) - 2*Pi) .lt. 1D-8) vtheta_bounds(2) = 0.0d0

            map_coords(4,:) = (/ grid_r(3), vtheta_bounds(2) /)
        end if
    else if (vr_prime .gt. (grid_r(ubound(grid_r,1)) - grid_r(2)/2)) then
        ! If the point is closer to lower vz bound.
        map_coords(1,:) = (/ grid_r(ubound(grid_r,1)), vtheta_bounds(1) /)
        map_coords(2,:) = (/ grid_r(ubound(grid_r,1)), vtheta_bounds(2) /)

        if (vtheta_bounds(2) .eq. 0.0d0) vtheta_bounds(2) = 2*Pi ! Set bound to 2*pi temporarily to ensure check works.
        if (vtheta_prime - vtheta_bounds(1) .lt. vtheta_bounds(2) - vtheta_prime) then
            ! If the point is closer to lower grid_theta bound.
            if (abs(vtheta_bounds(2) - 2*Pi) .lt. 1D-8) vtheta_bounds(2) = 0.0d0

            map_coords(3,:) = (/ grid_r(ubound(grid_r,1) - 1), vtheta_bounds(1) /)

            if (vtheta_bounds(1) .eq. 0.0d0) then
                map_coords(4,:) = (/ grid_r(ubound(grid_r,1)), grid_theta(ubound(grid_theta, 1)) /)
            else
                map_coords(4,:) = (/ grid_r(ubound(grid_r,1)), grid_theta(find_loc(grid_theta, vtheta_bounds(1)) - 1) /)
            end if
        else
            ! If the point is closer to higher grid_theta bound.
            if (abs(vtheta_bounds(2) - 2*Pi) .lt. 1D-8) vtheta_bounds(2) = 0.0d0

            map_coords(3,:) = (/ grid_r(ubound(grid_r,1) - 1), vtheta_bounds(2) /)

            if (vtheta_bounds(2) .eq. grid_theta(ubound(grid_theta,1))) then
                map_coords(4,:) = (/ grid_r(ubound(grid_r,1)), 0.0d0 /)
            else
                map_coords(4,:) = (/ grid_r(ubound(grid_r,1)), grid_theta(find_loc(grid_theta, vtheta_bounds(2)) + 1) /)
            end if
        end if
    else
        vr_bounds = find_v_bounds(vr_prime, grid_r)

        ! If the point is closer to lower vz bound.
        if (vr_prime - vr_bounds(1) .lt. vr_bounds(2) - vr_prime) then
            ! If the point is closer to the lower vr bound.
            map_coords(1,:) = (/ vr_bounds(1), vtheta_bounds(1) /)
            map_coords(2,:) = (/ vr_bounds(1), vtheta_bounds(2) /)

            if (vtheta_bounds(2) .eq. 0.0d0) vtheta_bounds(2) = 2*Pi ! Set bound to 2*pi temporarily to ensure check works.
            if (vtheta_prime - vtheta_bounds(1) .lt. vtheta_bounds(2) - vtheta_prime) then
                ! If the point is closer to the lower vtheta bound.
                if (abs(vtheta_bounds(2) - 2*Pi) .lt. 1D-8) vtheta_bounds(2) = 0.0d0

                map_coords(3,:) = (/ vr_bounds(2), vtheta_bounds(1) /)
                map_coords(4,:) = (/ grid_r(find_loc(grid_r, vr_bounds(1)) - 1), vtheta_bounds(1) /)
            else
                ! If the point is closer to the higher vtheta bound.
                if (abs(vtheta_bounds(2) - 2*Pi) .lt. 1D-8) vtheta_bounds(2) = 0.0d0

                map_coords(3,:) = (/ vr_bounds(2), vtheta_bounds(2) /)
                map_coords(4,:) = (/ grid_r(find_loc(grid_r, vr_bounds(1)) - 1), vtheta_bounds(2) /)
            end if
        else
            ! If the point is closer to the higher vr bound.
            map_coords(1,:) = (/ vr_bounds(2), vtheta_bounds(1) /)
            map_coords(2,:) = (/ vr_bounds(2), vtheta_bounds(2) /)

            if (vtheta_bounds(2) .eq. 0.0d0) vtheta_bounds(2) = 2*Pi ! Set bound to 2*pi temporarily to ensure check works.
            if (vtheta_prime - vtheta_bounds(1) .lt. vtheta_bounds(2) - vtheta_prime) then
                if (abs(vtheta_bounds(2) - 2*Pi) .lt. 1D-8) vtheta_bounds(2) = 0.0d0

                map_coords(3,:) = (/ vr_bounds(1), vtheta_bounds(1) /)
                map_coords(4,:) = (/ grid_r(find_loc(grid_r, vr_bounds(2)) + 1), vtheta_bounds(1) /)
            else
                if (abs(vtheta_bounds(2) - 2*Pi) .lt. 1D-8) vtheta_bounds(2) = 0.0d0

                map_coords(3,:) = (/ vr_bounds(1), vtheta_bounds(2) /)
                map_coords(4,:) = (/ grid_r(find_loc(grid_r, vr_bounds(2)) + 1), vtheta_bounds(2) /)
            end if
        end if
    end if
end subroutine find_points

! Given coordinates to map mass back to, calculate mass and add back to vdf.
subroutine replenish(map_coords, vdf, delta_m, grid_r, grid_theta, vr_prime, vtheta_prime, &
                     mass_sum, negative_count, positive_count)
    implicit none
    double precision, intent(in) :: map_coords(4,2)
    double precision, allocatable, intent(inout) :: vdf(:,:)
    double precision, allocatable, intent(in) :: grid_r(:), grid_theta(:)
    double precision, intent(in) :: vr_prime, vtheta_prime
    double precision, intent(in) :: delta_m

    double precision :: mapping(4,4) ! matrix A
    double precision :: b(4) ! matrix B
    double precision :: pivot(4) = 0.0
    integer :: rc ! return code
    integer :: i
    integer :: vr_loc, vtheta_loc ! index of v,grid_theta in vdf
    double precision :: total_p_x, total_p_y, total_e
    double precision, intent(inout) :: mass_sum
    integer, intent(inout) :: negative_count, positive_count

    ! Create mapping matrix.
    ! Mass.
    mapping(1,:) = 1.0d0
    b(1) = 1.0d0
    ! x-momentum.
    do i = 1,4
        mapping(2,i) = map_coords(i,1)*cos(map_coords(i,2))
    end do
    b(2) = vr_prime*cos(vtheta_prime)
    ! y-momentum.
    do i = 1,4
        mapping(3,i) = map_coords(i,1)*sin(map_coords(i,2))
    end do
    b(3) = vr_prime*sin(vtheta_prime)
    ! Energy.
    do i = 1,4
        mapping(4,i) = map_coords(i,1)**2
    end do
    b(4) = vr_prime**2
    b = delta_m * b

    ! Sanitize mapping.
    where(abs(mapping) .lt. 1D-8)
        mapping = 0.0
    endwhere

    ! LAPACK matrix solver.
    call dgesv(4, 1, mapping, 4, pivot, b, 4, rc)
    if (rc .ne. 0) then
        print *, "error on 1: ", rc
        do i = 1,4
            print *, map_coords(i,1), map_coords(i,2)
        end do
        print *, vr_prime, vtheta_prime
        print *, ""
        return
    end if

    ! Check remapping conserves momentum.
    total_p_x = 0.0d0
    total_p_y = 0.0d0
    total_e = 0.0d0

    do i = 1,5
        total_p_x = total_p_x + b(i) * map_coords(i,1)*cos(map_coords(i,2))
        total_p_y = total_p_y + b(i) * map_coords(i,1)*sin(map_coords(i,2))
        total_e = total_e + b(i) * map_coords(i,1)**2
    end do

    if (abs(total_p_x) - abs(delta_m*vr_prime*cos(vtheta_prime)) .gt. 1D-12) then 
        print *, "x_momentum not conserved"
        do i = 1,5
            print *, map_coords(i,:)
        end do
    end if
    if (abs(total_p_y) - abs(delta_m*vr_prime*sin(vtheta_prime)) .gt. 1D-12) then
        print *, "y_momentum not conserved"
        do i = 1,5
            print *, map_coords(i,:)
        end do
    end if
    if (abs(total_e) - abs(delta_m*vr_prime**2) .gt. 1D-12) then
        print *, "energy not conserved"
        print *, total_e, delta_m*vr_prime**2
    end if

    ! Add mass back to vdf.
    do i = 1,4
        vr_loc = find_loc(grid_r, map_coords(i,1))
        vtheta_loc = find_loc(grid_theta, map_coords(i,2))
        vdf(vr_loc, vtheta_loc) = vdf(vr_loc, vtheta_loc) + b(i)

        if (vr_loc .eq. 1) then
            mass_sum = mass_sum + b(i)
            if (b(i) .lt. 0.0) negative_count = negative_count + 1
            if (b(i) .gt. 0.0) positive_count = positive_count + 1
        end if
    end do
end subroutine replenish

! ------------------ COLLIDE SUBROUTINE HELPER FUNCTIONS ------------------ !
! Find grid angles that bound theta_v
pure function find_theta_bounds(theta_v, grid_theta) result(result)
    implicit none
    double precision, intent(in) :: theta_v
    double precision, allocatable, intent(in) :: grid_theta(:)
    double precision :: result(2)
    integer :: i

    result = 0.0

    if (theta_v .gt. grid_theta(ubound(grid_theta,dim=1))) then 
        result(1) = grid_theta(ubound(grid_theta,dim=1))
        result(2) = 0.0
    else
        do i = 1,size(grid_theta)-1
            if (theta_v .lt. grid_theta(i+1)) then
                result(1) = grid_theta(i)
                result(2) = grid_theta(i+1)
                exit
            end if
        end do
    end if
end function find_theta_bounds

! Find grid velocities that bound vr
pure function find_v_bounds(vr, grid_r) result(result)
    implicit none
    double precision, intent(in) :: vr
    double precision, allocatable, intent(in) :: grid_r(:)
    double precision :: result(2)
    integer :: i

    result = 0.0

    do i = 1,ubound(grid_r,1)-1
        if (vr .lt. grid_r(i+1) .or. abs(vr - grid_r(i+1)) .lt. 1D-12) then
            result(1) = grid_r(i)
            result(2) = grid_r(i+1)
            exit
        end if
    end do
end function find_v_bounds

! Better version of inbuilt findloc designed to compare floats.
pure function find_loc(arr, target) result(res)
    implicit none
    double precision, allocatable, intent(in) :: arr(:)
    double precision, intent(in) :: target
    integer :: res
    integer :: i

    do i = 1,size(arr)
        if (abs(arr(i) - target) .lt. 1D-12) res = i
    end do

    if (res .eq. 0) then
        res = 1
    end if

end function find_loc

! Binary search through cumulative distribution function to find pre-collision velocities.
pure function binary_search(cdf, target) result(l)
    implicit none
    double precision, allocatable, intent(in) :: cdf(:)
    double precision, intent(in) :: target
    integer :: low, high, mid
    integer :: l
    logical :: done

    done = .false.
    low = 1
    high = ubound(cdf,1)

    if (target .lt. cdf(1)) then
        l = 1
        done = .true.
    end if

    do while(.not. done)
        mid = (low + high)/2

        if (high - low == 1) then
            l = high
            done = .true.
        else if (target .gt. cdf(mid)) then
            low = mid
        else if (target .lt. cdf(mid)) then
            high = mid
        end if
    end do
end function binary_search

! ------------------ VDF HELPER FUNCTIONS ------------------ !
pure function calc_x_momentum(grid_r, grid_theta, vdf) result(x_momentum)
    implicit none
    double precision, allocatable, intent(in) :: grid_r(:)
    double precision, allocatable, intent(in) :: grid_theta(:)
    double precision, allocatable, intent(in) :: vdf(:,:)
    double precision :: x_momentum
    integer :: i, j

    x_momentum = 0.0
    do i = 1,size(vdf,1)
        do j = 1,size(vdf,2)
            x_momentum = x_momentum + vdf(i, j) * grid_r(i) * cos(grid_theta(j))
        end do
    end do
end function calc_x_momentum

pure function calc_y_momentum(grid_r, grid_theta, vdf) result(y_momentum)
    implicit none
    double precision, allocatable, intent(in) :: grid_r(:)
    double precision, allocatable, intent(in) :: grid_theta(:)
    double precision, allocatable, intent(in) :: vdf(:,:)
    double precision :: y_momentum
    integer :: i, j

    y_momentum = 0.0
    do i = 1,size(vdf,1)
        do j = 1,size(vdf,2)
            y_momentum = y_momentum + vdf(i, j) * grid_r(i) * sin(grid_theta(j))
        end do
    end do
end function calc_y_momentum

function calc_energy(grid_r, grid_theta, vdf) result(energy)
    implicit none
    double precision, allocatable, intent(in) :: grid_r(:)
    double precision, allocatable, intent(in) :: grid_theta(:)
    double precision, allocatable, intent(in) :: vdf(:,:)
    double precision :: energy
    integer :: i, j
    double precision :: vector(2)

    energy = 0.0
    do i = 1,size(vdf,1)
        do j = 1,size(vdf,2)
            vector = (/ grid_r(i)*cos(grid_theta(j)), grid_r(i)*sin(grid_theta(j)) /)
            energy = energy + vdf(i, j) * dot_product(vector, vector)
        end do
    end do
end function calc_energy

function calc_entropy(vdf) result(entropy)
    implicit none
    double precision, allocatable, intent(in) :: vdf(:,:)
    double precision :: entropy
    integer :: i, j

    entropy = 0.0
    do i = 1,size(vdf,1)
        do j = 1,size(vdf,2)
            if (vdf(i, j) .ne. 0.0d0) entropy = entropy + abs(vdf(i, j)) * log(abs(vdf(i, j)))
        end do
    end do
    entropy = entropy * -(1)/1.0d0

end function calc_entropy

end module polar_collide
