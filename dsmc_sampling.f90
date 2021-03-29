program dsmc_sampling
    use polar_collide
    implicit none

    double precision, parameter :: vr_min = 0.0d0
    double precision, parameter :: vr_max = 3.0d0
    double precision, parameter :: vz_min = -3.0d0
    double precision, parameter :: vz_max = 3.0d0
    integer, parameter :: n_r = 20 ! number of radial velocity grid points
    integer, parameter :: n_theta = 40 ! number of grid_theta grid points
    integer, parameter :: n_z = 40

    double precision, allocatable :: grid_r(:) ! radial velocities
    double precision, allocatable :: grid_theta(:) ! values from 0 to 360 - delta_theta converted to radians
    double precision, allocatable :: grid_z(:)
    double precision, allocatable :: vdf(:,:,:) ! velocity distribution function
    double precision :: particle_bin(20) ! count how many times particles fall into certain areas
    double precision :: v_bounds(2)

    double precision :: random_vx, random_vy, random_vz, random_vr, random_vtheta
    double precision :: delta_m, dr, dtheta, dz
    double precision :: map_coords(5,3)
    integer :: i, j, idx
    double precision :: x
    integer :: y, z

    allocate(grid_r(n_r))
    allocate(grid_theta(n_theta))
    allocate(grid_z(n_z))
    allocate(vdf(n_r, n_theta, n_z))

    grid_r(1) = vr_min
    do i = 2,n_r
        grid_r(i) = grid_r(i-1) + (vr_max - vr_min)/(n_r - 1) 
    end do
    dr = grid_r(2) - grid_r(1)

    grid_theta(1) = 0.0d0
    do i = 2,n_theta
        grid_theta(i) = grid_theta(i-1) + (2*Pi)/n_theta
    end do
    dtheta = grid_theta(2) - grid_theta(1)

    grid_z(1) = vz_min
    do i = 2,n_z
        grid_z(i) = grid_z(i-1) + (vz_max - vz_min)/(n_z - 1)
    end do
    dz = grid_z(2) - grid_z(1)

    vdf = 0.0d0
    delta_m = 1.5e-7

    do i = 1,10000000
        call random_init(.false., .true.)
        call random_number(random_vx)
        call random_number(random_vy)
        call random_number(random_vz)
        random_vx = random_vx * 6.4d0 - 3.2d0
        random_vy = random_vy * 6.4d0 - 3.2d0
        random_vz = random_vz * 6.4d0 - 3.2d0

        random_vr = sqrt(random_vx**2 + random_vy**2)
        random_vtheta = atan2(random_vy, random_vx)
        if (random_vtheta .lt. 0) random_vtheta = random_vtheta + 2*Pi
        if (abs(random_vr) .gt. vz_max .and. random_vr .gt. vr_max) cycle

        if (random_vz .gt. grid_z(20) - dz/2 .and. random_vz .lt. grid_z(20) + dz/2) then 
            if (random_vr .lt. dr/2) then
                particle_bin(1) = particle_bin(1) + 1
            else if (random_vr .gt. vr_max) then
                particle_bin(20) = particle_bin(20) + 1
            else
                v_bounds = find_v_bounds(random_vr, grid_r)
                if (random_vr - v_bounds(1) .lt. v_bounds(2) - random_vr) then
                    idx = find_loc(grid_r, v_bounds(1))
                    particle_bin(idx) = particle_bin(idx) + 1
                else
                    idx = find_loc(grid_r, v_bounds(2))
                    particle_bin(idx) = particle_bin(idx) + 1
                end if
            end if
        end if

        call find_points(random_vr, random_vtheta, random_vz, grid_r, grid_theta, grid_z, map_coords)
        ! do j = 1,5
        !     print *, map_coords(j,:)
        ! end do
        ! print *, random_vr, random_vtheta, random_vz
        ! print *, ""
        ! if (random_vr .lt. (grid_r(2) + grid_r(3))/2.0d0) then
        !     do j = 1,4
        !         print *, map_coords(j,:)
        !     end do
        ! end if
        call replenish(map_coords, vdf, delta_m, grid_r, grid_theta, grid_z, random_vr, random_vtheta, random_vz, &
        x, y, z)
        
    end do

    open(unit=20, file="dsmc.txt")
    do i = 1,n_r
        write(20,*) vdf(i,:,5)
    end do
    close(20)

    open(unit=21, file="dsmc.dat", access="stream")
    write(21) vdf
    close(21)

    open(unit=22, file="particle_bins.dat", access="stream")
    write(22) particle_bin
    close(22)

    print *, sum(vdf)
    print *, particle_bin

    deallocate(grid_r)
    deallocate(grid_theta)
    deallocate(vdf)

end program dsmc_sampling

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
