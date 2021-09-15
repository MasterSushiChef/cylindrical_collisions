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
    double precision, parameter :: vr_max = 4.5d0
    integer, parameter :: n_r = 9  ! number of radial velocity grid points
    integer, parameter :: n_theta = 28 ! number of theta grid points
    integer, parameter :: n_t = 200 ! number of timesteps
    double precision, parameter :: m_hat = 1
    double precision, parameter :: t_hat = 0.1d0
    double precision, parameter :: ndens_hat = 1
    double precision, parameter :: temp_hat = 1.0
    double precision, parameter :: kn = 1.0d0
    
    double precision, parameter :: k0 = 0.6d0
    double precision, parameter :: crms = 5D-3

    ! select method used to deplete: 0 - N^2, 1 - Full Monte Carlo, 2 - Monte Carlo/N^2 Hybrid
    integer, parameter :: method = 1
    character(len = 8), parameter :: fmt = '(I6.6)'
    character(len = 6) :: x1
    character(len = 15) :: file_name

    integer :: i, j, t, x, y, vr, vtheta ! loop variables

    double precision, allocatable :: grid_r(:) ! radial velocities
    double precision, allocatable :: grid_theta(:) ! values from 0 to 360 - delta_theta converted to radians
    double precision, allocatable :: vdf(:,:) ! velocity distribution function
    double precision, allocatable :: cdf(:) ! cumulative distribution function
    ! double precision, allocatable :: equal_area_remapping(:) ! points to separate equal area remappings

    double precision, allocatable :: mass_collector(:,:)
    double precision :: initial_zero_point

    double precision :: map_coords(5,2) ! coordinates on velocity grid to map post collision mass back onto
    double precision :: dr, dtheta
    double precision :: nc ! number of collisions
    double precision :: n_hat_neg

    double precision :: vr1, vr2, vtheta1, vtheta2! pre collision velocities
    double precision :: vr1_prime, vr2_prime, vtheta1_prime, vtheta2_prime ! post collision velocities
    double precision :: delta_m1, delta_m2 ! colliding mass
    double precision :: mass1, x_momentum1, y_momentum1, energy1
    double precision :: mass2, x_momentum2, y_momentum2, energy2
    double precision :: entropy(n_t+1), moment(n_t+1), neg_mass(n_t+1), zero_point(n_t+1)
    integer :: cutoff

    allocate(grid_r(n_r))
    allocate(grid_theta(n_theta))
    ! allocate(equal_area_remapping(n_r-1))
    allocate(vdf(n_r, n_theta))
    allocate(mass_collector(n_r, n_theta))
    allocate(cdf(n_r * n_theta - (n_theta - 1)))

    ! Initialize matrices to collect statistics.
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
    vdf(1,1) = exp(0.0d0) * (m_hat/temp_hat) * Pi * (dr/2)**2
    do vr = 2,size(grid_r)
        do vtheta = 1,size(grid_theta)
            vdf(vr, vtheta) = exp(-((grid_r(vr)*cos(grid_theta(vtheta)))**2 + &
                (grid_r(vr)*sin(grid_theta(vtheta)))**2) * (m_hat/temp_hat))
            vdf(vr, vtheta) = vdf(vr, vtheta) * grid_r(vr) * dr * dtheta
        end do
    end do
    vdf = vdf * ndens_hat * (m_hat/(Pi*temp_hat))**1.0
    vdf = vdf/sum(vdf)
    initial_zero_point = vdf(1,1)

    if (method .eq. 2) then 
        cutoff = n_r/2 ! last vr to do Monte Carlo, switch to N^2 afterwards
    else
        cutoff = 0
    end if

    ! vdf(2,1) = 0.5d0
    ! vdf(2,21) = 0.5d0

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
    ! vdf = vdf/sum(vdf)

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
    moment(1) = calc_moment(vdf, grid_r, 8)
    neg_mass(1) = abs(sum(vdf, mask=vdf .lt. 0.0d0))/sum(vdf)
    zero_point(1) = vdf(1, 1)
    print *, "Initial mass: ", mass1
    print *, "Initial x-momentum: ", x_momentum1
    print *, "Initial y-momentum: ", y_momentum1
    print *, "Initial energy: ", energy1
    print *, ""

    do t = 1,n_t
        if (method .eq. 0 .or. method .eq. 2) then
            do y = 1,n_theta
                do x = cutoff+1,n_r
                    do j = 1,n_theta
                        do i = cutoff+1,n_r
                            ! Calculate delta_m. Psuedo-Maxwell molecules.
                            if (x .eq. i .and. y .eq. j) cycle
                            if (vdf(i, j) .eq. 0.0d0 .or. vdf(x, y) .eq. 0.0d0) cycle
                            delta_m1 = 0.5 * t_hat * vdf(x, y) * vdf(i, j)
                            delta_m2 = delta_m1

                            ! Pre-collision velocities.
                            vr1 = grid_r(x)
                            vtheta1 = grid_theta(y)
                            vr2 = grid_r(i)
                            vtheta2 = grid_theta(j)

                            ! Collide two particles.
                            call collide(vr1, vr2, vtheta1, vtheta2, vr1_prime, vr2_prime, vtheta1_prime, vtheta2_prime)

                            ! Deplete from vdf.
                            call deplete(vdf, grid_r, grid_theta, vr1, vr2, vtheta1, vtheta2, delta_m1, delta_m2)

                            ! Find points to map mass back to and add it to vdf for both points.
                            call find_points(vr1_prime, vtheta1_prime, grid_r, grid_theta, map_coords)
                            call replenish(map_coords, vdf, delta_m1, grid_r, grid_theta, vr1_prime, vtheta1_prime)
                            call find_points(vr2_prime, vtheta2_prime, grid_r, grid_theta, map_coords)
                            call replenish(map_coords, vdf, delta_m2, grid_r, grid_theta, vr2_prime, vtheta2_prime)
                        end do
                    end do
                end do
            end do

            print *, t

            entropy(t+1) = calc_entropy(vdf)
            moment(t+1) = calc_moment(vdf, grid_r, 8)
            neg_mass(t+1) = abs(sum(vdf, mask=vdf .lt. 0.0d0))/sum(vdf)

            ! Write out vdf data for post processing.
            write (x1, fmt) t
            file_name = "vdf_" // trim(x1) // ".dat"
            open(unit=20, file=file_name, access="stream")
            write(20) vdf
            close(20)
        end if
        if (method .eq. 1 .or. method .eq. 2) then
            ! Build cdf.
            cdf = 0.0d0
            call build_cdf(n_r, n_theta, vdf, cdf)
            nc = nint((t_hat * temp_hat**(2.0/3.0d0))/(kn * crms**2 * 1.0d0 * sum(grid_r * dr * dtheta)/n_r)) ! beta^3 avg

            ! Calculate delta_m. Using psuedo-Maxwell molecules.
            n_hat_neg = sum(vdf, mask=vdf .lt. 0.0d0)
            delta_m1 = (t_hat * (sum(vdf) - 2*n_hat_neg)**2)/(2.0d0 * kn * nc)
            delta_m2 = delta_m1
            print *, delta_m1, nc, t

            do j = 1,int(nc)
                ! Calculate pre-collision velocities.
                call precollision(grid_r, grid_theta, n_theta, cdf, vr1, vtheta1)
                call precollision(grid_r, grid_theta, n_theta, cdf, vr2, vtheta2)

                ! If we select the same velocity, cycle.
                if ((vr1 .eq. vr2) .and. (vtheta1 .eq. vtheta2)) cycle

                ! Collide two particles.
                call collide(vr1, vr2, vtheta1, vtheta2, vr1_prime, vr2_prime, vtheta1_prime, vtheta2_prime)
                
                ! Multiply delta_m by signs to determine depletion direction and deplete. (full Monte Carlo).
                call deplete(vdf, grid_r, grid_theta, vr1, vr2, vtheta1, vtheta2, delta_m1, delta_m2)

                ! Find points to map mass back to and add it to vdf for both points.
                call find_points(vr1_prime, vtheta1_prime, grid_r, grid_theta, map_coords)
                call replenish(map_coords, vdf, delta_m1, grid_r, grid_theta, vr1_prime, vtheta1_prime)
                call find_points(vr2_prime, vtheta2_prime, grid_r, grid_theta, map_coords)
                call replenish(map_coords, vdf, delta_m2, grid_r, grid_theta, vr2_prime, vtheta2_prime)
            end do

            entropy(t+1) = calc_entropy(vdf)
            moment(t+1) = calc_moment(vdf, grid_r, 8)
            neg_mass(t+1)  = abs(sum(vdf, mask=vdf .lt. 0.0d0))/sum(vdf)
            zero_point(t+1) = vdf(1, 1)

            ! Write out vdf data for post processing.
            write (x1, fmt) t
            file_name = "vdf_" // trim(x1) // ".dat"
            open(unit=20, file=file_name, access="stream")
            write(20) vdf
            close(20)
        end if
        if (method .gt. 2) then
            print *, "Unknown method selected."
            stop
        end if
    end do

    ! Write data out for post processing in Python.
    open(unit=21, file="vdf.txt", action="write", status="old")
    do i = 1,n_r
        write(21,*) vdf(i,:)
    end do

    open(unit=22, file="entropy_theta50_6.dat", access="stream")
    write(22) entropy
    close(22)

    open(unit=25, file="M8_maxwellian_working.dat", access="stream")
    write(25) moment
    close(25)

    open(unit=23, file="negative_mass_theta50_6.dat", access="stream")
    write(23) neg_mass
    close(23)

    open(unit=26, file="zero_point_theta50_6.dat", access="stream")
    write(26) zero_point
    close(26)

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
    print *, "Zero point initial: ", initial_zero_point
    print *, "Zero point after: ", vdf(1, 1)
    print *, "Zero point percent change: ", (vdf(1,1) - initial_zero_point)/initial_zero_point * 100
    print *, "Negative mass in vdf: ", sum(vdf, mask=vdf .lt. 0.0d0)

    deallocate(vdf)
    deallocate(grid_r)
    deallocate(grid_theta)
    deallocate(cdf)

end program collision
