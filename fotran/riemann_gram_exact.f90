! Exact Gram data on the critical line, from Proposition 4.1 (check of riemann_dn):
!
!     <G_j, G_k>_{1/2} = (1/(jk)) sum_{m>=1} (m mod j)(m mod k) / (m(m+1)),
!     <E, G_k>_{1/2}   = -(1/k)   sum_{m>=1} (m mod k) / (m(m+1)) = -log(k)/k,
!
! where the series are summed exactly over one period L = lcm(j,k) by means of
!     sum_{q>=0} 1/((qL+r)(qL+r+1)) = (psi((r+1)/L) - psi(r/L)) / L,   1 <= r < L.
! No zeta function and no quadrature are involved.  As in riemann_dn, the Cholesky factor
! of the Gram matrix gives d_n(1/2)^2 for every n <= nmax, and s_n(1/2) is computed from
! the same Gram matrix.  The cost is O(lcm(j,k)) per entry, about nmax^4/10 in total.
!
! Usage:  ./riemann_gram_exact [file.nml]     (namelist &params, nmax, nthreads, out_file)

program riemann_gram_exact
    use riemann_lib
    use, intrinsic :: iso_fortran_env, only: int64
    !$ use omp_lib, only: omp_set_num_threads
    implicit none

    external :: dpotrf, dtrsv

    integer :: nmax, nthreads
    character(len=256) :: out_file, nml_file
    namelist /params/ nmax, nthreads, out_file

    real(dp), allocatable :: A(:,:), b(:), y_(:), d2(:), sg(:)
    integer,  allocatable :: mu(:), spf(:)
    integer  :: nk, j, k, n, n1, info, u, nok
    integer(int64) :: clock0, clock1, rate
    real(dp) :: acc, ttotal, berr
    character(len=:), allocatable :: line

    nmax     = 640
    nthreads = 20
    out_file = "dn_exact.txt"

    if (command_argument_count() >= 1) then
        call get_command_argument(1, nml_file)
        open(newunit=u, file=trim(nml_file), status="old", action="read")
        read(u, nml=params)
        close(u)
    end if
    !$ if (nthreads > 0) call omp_set_num_threads(nthreads)
    call system_clock(clock0, rate)

    ! ---- exact Gram matrix (upper triangle) and b_k = -log(k)/k ----
    nk = nmax - 1
    allocate(A(nk, nk), b(nk))
    A = 0.0_dp
    !$omp parallel do schedule(dynamic, 1) private(j)
    do k = nmax, 2, -1
        do j = 2, k
            A(j - 1, k - 1) = period_sum(j, k) / (real(j, dp) * real(k, dp))
        end do
    end do
    !$omp end parallel do

    ! check of the same period sums on b: sum_m (m mod k)/(m(m+1)) = log k
    berr = 0.0_dp
    do k = 2, nmax
        b(k - 1) = -log(real(k, dp)) / real(k, dp)
        berr = max(berr, abs(period_sum1(k) / log(real(k, dp)) - 1.0_dp))
    end do

    ! ---- s_n(1/2) = mu^T A mu - 2 mu^T b + 1 from the same Gram matrix ----
    allocate(mu(nmax), spf(nmax), sg(nmax))
    call sieve(nmax, mu, spf)
    sg(1) = 1.0_dp
    do n = 2, nmax
        sg(n) = sg(n - 1)
        if (mu(n) /= 0) then
            acc = real(mu(n), dp) * A(n - 1, n - 1)
            do k = 2, n - 1
                if (mu(k) /= 0) acc = acc + 2.0_dp * real(mu(k), dp) * A(k - 1, n - 1)
            end do
            sg(n) = sg(n) + real(mu(n), dp) * acc - 2.0_dp * real(mu(n), dp) * b(n - 1)
        end if
    end do

    ! ---- Cholesky A = U^T U, y = U^{-T} b, d_n^2 = 1 - sum_{i<n} y_i^2 ----
    call dpotrf("U", nk, A, nk, info)
    nok = nk
    if (info > 0) nok = info - 1
    allocate(y_(nok), d2(nok + 1))
    y_ = b(1:nok)
    call dtrsv("U", "T", "N", nok, A, nk, y_, 1)
    d2(1) = 1.0_dp
    do n1 = 1, nok
        d2(n1 + 1) = d2(n1) - y_(n1)**2
    end do

    call system_clock(clock1)
    ttotal = real(clock1 - clock0, dp) / real(rate, dp)

    open(newunit=u, file=trim(out_file), status="replace", action="write")
    write(u, '(a)') "# exact d_n(1/2)^2 and s_n(1/2) from the Gram matrix of Proposition 4.1"
    write(u, '(a,i0,a,i0,a,es9.2,a,f10.2,a)') "# nmax = ", nmax, ", stable up to n = ", nok + 1, &
        ", max |sum (m mod k)/(m(m+1)) / log k - 1| = ", berr, ", wall time = ", ttotal, " s"
    write(u, '(a)') "# n d_n^2 s_n"
    do n = 2, nok + 1
        line = inum(n) // " " // num(d2(n), '(es23.15)') // " " // num(sg(n), '(es23.15)')
        write(u, '(a)') line
    end do
    close(u)

    print '(a,i0,a,es9.2,a,f10.2,a,a)', "done: n <= ", nok + 1, ", b check ", berr, &
        ", wall time ", ttotal, " s  ->  ", trim(out_file)

contains

    ! sum_{m>=1} (m mod j)(m mod k) / (m(m+1)), exactly over one period L = lcm(j,k)
    real(dp) function period_sum(j, k) result(ssum)
        integer, intent(in) :: j, k

        integer  :: L, r, rj, rk, g, a, c
        real(dp) :: psi0, psi1, rL

        a = j
        c = k
        do while (c /= 0)
            g = mod(a, c)
            a = c
            c = g
        end do
        L  = (j / a) * k
        rL = real(L, dp)
        ssum = 0.0_dp
        psi0 = digamma(1.0_dp / rL)
        rj = 0
        rk = 0
        do r = 1, L - 1
            rj = rj + 1
            if (rj == j) rj = 0
            rk = rk + 1
            if (rk == k) rk = 0
            psi1 = digamma(real(r + 1, dp) / rL)
            ssum = ssum + real(rj, dp) * real(rk, dp) * (psi1 - psi0)
            psi0 = psi1
        end do
        ssum = ssum / rL
    end function period_sum

    ! sum_{m>=1} (m mod k) / (m(m+1)), exactly over one period
    real(dp) function period_sum1(k) result(ssum)
        integer, intent(in) :: k

        integer  :: r
        real(dp) :: rk

        rk = real(k, dp)
        ssum = 0.0_dp
        do r = 1, k - 1
            ssum = ssum + real(r, dp) * (digamma(real(r + 1, dp) / rk) - digamma(real(r, dp) / rk))
        end do
        ssum = ssum / rk
    end function period_sum1

end program riemann_gram_exact
