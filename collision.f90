! ------------------ MAIN PROGRAM ------------------ !
! Particle collisions in 2D using a polar coordinate system.
! zero velocity point.
! Problem is definitely near center points, edge stability is greatly affected by n_theta.
! Try finding proper split in middle using n_r = 20 and n_theta = 24.
! Experiement with other remappings.
program collision
    use polar_collide
    implicit none

    ! Declare variables.
    double precision, parameter :: vr_min = 0.0d0
    double precision, parameter :: vr_max = 3.0d0
    integer, parameter :: n_r = 20 ! number of radial velocity grid points
    integer, parameter :: n_theta = 68 ! number of theta grid points
    integer, parameter :: n_t = 250 ! number of timesteps
    integer, parameter :: m_hat = 1
    double precision, parameter :: t_hat = 0.1d0
    integer, parameter :: ndens_hat = 1
    integer, parameter :: temp_hat = 1
    double precision, parameter :: kn = 1.0d0
    
    double precision, parameter :: k0 = 0.6d0
    double precision, parameter :: crms = 5D-3

    character(len = 8), parameter :: fmt = '(I6.6)'
    character(len = 6) :: x1
    character(len = 15) :: file_name

    integer :: i, j, k, vr, vtheta ! loop variables

    double precision, allocatable :: grid_r(:) ! radial velocities
    double precision, allocatable :: grid_theta(:) ! values from 0 to 360 - delta_theta converted to radians
    double precision, allocatable :: vdf(:,:) ! velocity distribution function
    double precision, allocatable :: cdf(:) ! cumulative distribution function
    ! double precision, allocatable :: equal_area_remapping(:) ! points to separate equal area remappings

    double precision, allocatable :: remap_count(:,:)
    double precision, allocatable :: depletion_count(:,:)
    double precision, allocatable :: mass_collector(:,:)

    double precision :: map_coords(4,2) ! coordinates on velocity grid to map post collision mass back onto
    double precision :: dr, dtheta
    double precision :: nc ! number of collisions
    integer :: sign_one, sign_two
    integer :: vr1_idx, vr2_idx, vtheta1_idx, vtheta2_idx
    integer :: vr_loc, vtheta_loc

    double precision :: vr1, vr2, vtheta1, vtheta2! pre collision velocities
    double precision :: vr1_prime, vr2_prime, vtheta1_prime, vtheta2_prime ! post collision velocities
    double precision :: delta_m1, delta_m2, delta_m ! colliding mass
    double precision :: mass1, x_momentum1, y_momentum1, energy1
    double precision :: mass2, x_momentum2, y_momentum2, energy2
    double precision :: entropy(n_t+1)
    double precision :: mass_sum, initial_zero_point
    integer :: negative_count, positive_count

    allocate(grid_r(n_r))
    allocate(grid_theta(n_theta))
    ! allocate(equal_area_remapping(n_r-1))
    allocate(vdf(n_r, n_theta))
    allocate(mass_collector(n_r, n_theta))
    allocate(remap_count(n_r, n_theta))
    allocate(depletion_count(n_r, n_theta))
    allocate(cdf(n_r * n_theta - (n_theta - 1)))

    ! Initialize matrices to collect statistics.
    remap_count = 0
    depletion_count = 0
    mass_collector = 0.0d0

    ! Initialize values of velocity grid (vr, vtheta, vz).
    grid_r(1) = vr_min
    do i = 2,n_r
        grid_r(i) = grid_r(i-1) + (vr_max - vr_min)/(n_r - 1) 
    end do

    grid_theta(1) = 0.0d0
    do i = 2,n_theta
        grid_theta(i) = grid_theta(i-1) + (2*Pi)/n_theta
    end do

    dr = grid_r(2) - grid_r(1)
    dtheta = grid_theta(2) - grid_theta(1)

    ! Calculate values to have equal area remappings.
    ! do i = 1,n_r-1
    !     equal_area_remapping(i) = sqrt(0.5 * grid_r(i+1)**2 + 0.5 * grid_r(i)**2)
    ! end do

    ! Build Maxwellian velocity distribution function.
    vdf = 0.0d0
    vdf(1,1) = exp(-(0.0d0)**2) * (m_hat/temp_hat) * Pi * (dr/2)**2
    ! Set vdf(1,1) as 0 point, have to multiply by circle area element, new finding points, make sure it can't remap to other thetas.
    do vr = 2,size(grid_r)
        do vtheta = 1,size(grid_theta)
            vdf(vr, vtheta) = exp(-((grid_r(vr))**2) * (m_hat/temp_hat))
            vdf(vr, vtheta) = vdf(vr, vtheta) * grid_r(vr) * dr * dtheta
        end do
    end do
    vdf = vdf * ndens_hat * (m_hat/(Pi*temp_hat))**1.0
    vdf = vdf/sum(vdf)
    initial_zero_point = vdf(1,1)

    ! vdf(20,20) = 1.0d0
    ! vdf(15,1) = 1.0d0

    ! Build BKW velocity distribution function.
    ! vdf(1,:) = 0.0d0
    ! vdf(1,1) = 1/((2*k0) * (Pi*k0**(1.5d0))) * (5*k0 - 3) * Pi * (dr/2)**2
    ! do vtheta = 1,size(grid_theta)
    !   do vr = 2,size(grid_r)
    !     vdf(vr, vtheta) = (1/((2*k0) * (Pi*k0**(1.5d0))) * (5*k0 - 3  + (2*(1-k0))/(k0) * (grid_r(vr)**2 + &
    !     (grid_r(vr)**2)) * exp(-((grid_r(vr)**2 + (grid_r(vr)**2)**2)/k0))))
    !     vdf(vr, vtheta) = vdf(vr, vtheta) * dr * dtheta * grid_r(vr) ! how accurate is r dr dtheta compare to real area?
    !   end do
    ! end do

    open(unit=20, file="vdf_000000.dat", access="stream")
    write(20) vdf
    close(20)

    ! call find_points(0.01d0, 5.9d0, vr_max, grid_r, grid_theta, map_coords, equal_area_remapping)
    ! do i = 1,4
    !   print *, map_coords(i,:)
    ! end do

    ! Calculate initial mass, momenta, energy in vdf.
    mass1 = sum(vdf)
    x_momentum1 = calc_x_momentum(grid_r, grid_theta, vdf)
    y_momentum1 = calc_y_momentum(grid_r, grid_theta, vdf)
    energy1 = calc_energy(grid_r, grid_theta, vdf)
    entropy(1) = calc_entropy(vdf)
    print *, "Initial mass: ", mass1
    print *, "Initial x-momentum: ", x_momentum1
    print *, "Initial y-momentum: ", y_momentum1
    print *, "Initial energy: ", energy1
    print *, ""

    do i = 1,n_t
        ! Build cdf.
        cdf = 0.0d0
        call build_cdf(n_r, n_theta, vdf, cdf)

        ! Calculate delta_m. Using psuedo-Maxwell molecules.
        nc = nint((t_hat * temp_hat**(2.0/3.0d0))/(kn * crms**2 * 1.0d0 * sum(grid_r * dr * dtheta)/n_r)) ! beta^3 avg
        delta_m1 = (t_hat * sum(vdf)**2)/(2.0d0 * kn * nc)
        delta_m2 = delta_m1
        delta_m = delta_m1
        print *, delta_m1, nc, i

        do j = 1,int(nc)
            ! Calculate pre-collision velocities.
            call precollision(grid_r, grid_theta, n_r, n_theta, cdf, vr1, vtheta1)
            call precollision(grid_r, grid_theta, n_r, n_theta, cdf, vr2, vtheta2)
           ! print *, ""

            ! If we select the same velocity, cycle.
            if ((vr1 .eq. vr2) .and. (vtheta1 .eq. vtheta2)) cycle

            ! Collide two particles.
            call collide(vr1, vr2, vtheta1, vtheta2, vr1_prime, vr2_prime, vtheta1_prime, vtheta2_prime)
            
            ! Multiply delta_m by signs to determine depletion direction and deplete. (full Monte Carlo).
            vr1_idx = find_loc(grid_r, vr1)
            vr2_idx = find_loc(grid_r, vr2)
            vtheta1_idx = find_loc(grid_theta, vtheta1)
            vtheta2_idx = find_loc(grid_theta, vtheta2)
            sign_one = int(sign(1.0d0, vdf(vr1_idx, vtheta1_idx)))
            sign_two = int(sign(1.0d0, vdf(vr2_idx, vtheta2_idx)))
            delta_m1 = abs(delta_m1) * sign_one * sign_two
            delta_m2 = abs(delta_m2) * sign_one * sign_two
            vdf(vr1_idx, vtheta1_idx) = vdf(vr1_idx, vtheta1_idx) - delta_m1
            vdf(vr2_idx, vtheta2_idx) = vdf(vr2_idx, vtheta2_idx) - delta_m2
            depletion_count(vr1_idx, vtheta1_idx) = depletion_count(vr1_idx, vtheta1_idx) + 1
            depletion_count(vr2_idx, vtheta2_idx) = depletion_count(vr2_idx, vtheta2_idx) + 1

            ! Find points to map mass back to and add it to vdf for both points.
            call find_points(vr1_prime, vtheta1_prime, grid_r, grid_theta, map_coords)
            call replenish(map_coords, vdf, delta_m1, grid_r, grid_theta, vr1_prime, vtheta1_prime, &
                           mass_sum, negative_count, positive_count)
            do k = 1,4
                vr_loc = find_loc(grid_r, map_coords(k,1))
                vtheta_loc = find_loc(grid_theta, map_coords(k,2))
                remap_count(vr_loc, vtheta_loc) = remap_count(vr_loc, vtheta_loc) + 1
            end do
            call find_points(vr2_prime, vtheta2_prime, grid_r, grid_theta, map_coords)
            call replenish(map_coords, vdf, delta_m2, grid_r, grid_theta, vr2_prime, vtheta2_prime, &
                           mass_sum, negative_count, positive_count)
            do k = 1,4
                vr_loc = find_loc(grid_r, map_coords(k,1))
                vtheta_loc = find_loc(grid_theta, map_coords(k,2))
                remap_count(vr_loc, vtheta_loc) = remap_count(vr_loc, vtheta_loc) + 1 
            end do
        end do

        write (x1, fmt) i
        file_name = "vdf_" // trim(x1) // ".dat"
        open(unit=20, file=file_name, access="stream")
        write(20) vdf
        close(20)

        entropy(i+1) = calc_entropy(vdf)
    end do

    ! Write data out for post processing in Python.
    open(unit=21, file="vdf.txt", action="write", status="old")
    do i = 1,n_r
        write(21,*) vdf(i,:)
    end do

    open(unit=24, file="depletion.txt", action="write", status="old")
    do i = 1,n_r
        write(24,*) depletion_count(i,:)
    end do

    open(unit=23, file="remap_count.txt", action="write", status="old")
    do i = 1,n_r
        write(23,*) remap_count(i,:)
    end do

    open(unit=22, file="entropy.dat", access="stream")
    write(22) entropy
    close(22)

    ! Check mass, momentum, energy are conserved.
    mass2 = sum(vdf)
    x_momentum2 = calc_x_momentum(grid_r, grid_theta, vdf)
    y_momentum2 = calc_y_momentum(grid_r, grid_theta, vdf)
    energy2 = calc_energy(grid_r, grid_theta, vdf)
    print *, ""
    print *, "Mass after collisions: ", mass2
    print *, "x-momentum after collisions: ", x_momentum2
    print *, "y-momentum after collisions: ", y_momentum2
    print *, "Energy after collisions: ", energy2
    print *, ""

    print *, "! ------------------ STATS ------------------ !"
    print *, "Mass percent error: ", (mass2 - mass1)/mass1 * 100
    print *, "x-momentum percent error: ", (x_momentum2 - x_momentum1)/x_momentum1 * 100
    print *, "y-momentum percent error: ", (y_momentum2 - y_momentum1)/y_momentum1 * 100
    print *, "Energy percent error: ", (energy2 - energy1)/energy1 * 100
    ! print *, "Net remapped mass on zero point: ", mass_sum
    print *, "Times chosen as internal point: ", positive_count
    print *, "Times chosen as external point: ", negative_count
    print *, "Zero point initial: ", initial_zero_point
    print *, "Zero point after: ", vdf(1, 1)

    deallocate(vdf)
    deallocate(grid_r)
    deallocate(grid_theta)
    deallocate(cdf)
    deallocate(depletion_count)
    deallocate(remap_count)

end program collision
