! Optimal error on the critical line (Subsection 4.3 of the paper) in Fortran:
!
!     d_n(1/2)^2 = min_c || E - sum_{k=2}^n c_k G_k ||^2_{1/2} = ||E||^2 - b^T A^{-1} b,
!
! with A = (<G_j, G_k>)_{j,k=2}^{nmax} and b = (<E, G_k>)_{k=2}^{nmax}.  Since E and the
! G_k are real on the real axis, <F, H> = (1/pi) int_0^inf Re[F(1/2+it) conj H(1/2+it)] dt.
! These integrals are evaluated on a fixed grid (uniform step h_fine on [0, T_fine),
! n_geo geometric nodes on [T_fine, T_max], trapezoidal weights), with zeta(1/2+it)
! from the Euler-Maclaurin formula (zeta_em in riemann_lib.f90).
!
! With A = U^T U (Cholesky) and y = U^{-T} b, the leading blocks give every n at once:
!     d_n^2 = ||E||^2 - sum_{i=1}^{n-1} y_i^2,   2 <= n <= nmax.
! The Moebius error s_n(1/2) = mu^T A mu - 2 mu^T b + ||E||^2 is also computed from the
! same Gram data, as a check against riemann_traj.
!
! Usage:  ./riemann_dn [file.nml]     (namelist &params, e.g. riemann_dn.nml)

program riemann_dn
    use riemann_lib
    use, intrinsic :: iso_fortran_env, only: int64
    !$ use omp_lib, only: omp_set_num_threads
    implicit none

    external :: dgemm, dpotrf, dtrsv

    ! ---- parameters (namelist &params); defaults reproduce the grid of the paper ----
    integer  :: nmax, n_geo, block, nthreads
    real(dp) :: T_fine, h_fine, T_max
    character(len=256) :: out_file, nml_file
    namelist /params/ nmax, T_fine, h_fine, T_max, n_geo, block, nthreads, out_file

    ! ---- work arrays ----
    real(dp), allocatable :: t(:), w(:), X(:,:), Y(:,:), A(:,:), b(:), y_(:), d2(:), sg(:)
    integer,  allocatable :: mu(:), spf(:)
    integer  :: nt, n1, nk, i, k, j0, j1, n, info, u, nfine, nok
    integer(int64) :: clock0, clock1, rate
    real(dp) :: E2, om, lk, acc, ttotal
    complex(dp) :: s, z, psi, ek, gk
    character(len=:), allocatable :: line

    nmax     = 5120
    T_fine   = 80.0_dp
    h_fine   = 0.01_dp
    T_max    = 3.0e5_dp
    n_geo    = 20000
    block    = 128
    nthreads = 20
    out_file = "dn_f90.txt"

    if (command_argument_count() >= 1) then
        call get_command_argument(1, nml_file)
        open(newunit=u, file=trim(nml_file), status="old", action="read")
        read(u, nml=params)
        close(u)
    end if
    !$ if (nthreads > 0) call omp_set_num_threads(nthreads)
    call system_clock(clock0, rate)

    ! ---- quadrature grid on [0, T_max] with trapezoidal weights ----
    nfine = nint(T_fine / h_fine)
    nt = nfine + n_geo
    allocate(t(nt), w(nt))
    t(1:nfine) = [(real(i - 1, dp) * h_fine, i = 1, nfine)]
    t(nfine+1:nt) = [(T_fine * (T_max / T_fine)**(real(i - 1, dp) / real(n_geo - 1, dp)), &
                      i = 1, n_geo)]
    w(1)  = 0.5_dp * (t(2) - t(1))
    w(nt) = 0.5_dp * (t(nt) - t(nt - 1))
    w(2:nt-1) = 0.5_dp * (t(3:nt) - t(1:nt-2))

    ! ---- sampled G_k: rows = nodes, columns = k = 2..nmax, scaled by sqrt(w/pi) ----
    nk = nmax - 1
    allocate(X(nt, nk), Y(nt, nk), b(nk))
    E2 = 0.0_dp
    b  = 0.0_dp
    !$omp parallel do schedule(dynamic, 16) private(s, z, psi, ek, gk, om, k, lk) &
    !$omp& reduction(+:E2, b)
    do i = 1, nt
        s   = cmplx(0.5_dp, t(i), dp)
        z   = zeta_em(s)
        psi = z / s
        ek  = 1.0_dp / s
        om  = w(i) / pi
        E2  = E2 + om * abs(ek)**2
        do k = 2, nmax
            lk = log(real(k, dp))
            gk = (exp(-s * lk) - 1.0_dp / real(k, dp)) * psi
            X(i, k - 1) = sqrt(om) * real(gk)
            Y(i, k - 1) = sqrt(om) * aimag(gk)
            b(k - 1) = b(k - 1) + om * (real(ek) * real(gk) + aimag(ek) * aimag(gk))
        end do
    end do
    !$omp end parallel do

    ! ---- Gram matrix A = X^T X + Y^T Y (upper triangle), by column blocks ----
    allocate(A(nk, nk))
    A = 0.0_dp
    !$omp parallel do schedule(dynamic, 1) private(j0, j1)
    do j0 = 1, nk, block
        j1 = min(nk, j0 + block - 1)
        call dgemm("T", "N", j1, j1 - j0 + 1, nt, 1.0_dp, X, nt, X(1, j0), nt, &
                   0.0_dp, A(1, j0), nk)
        call dgemm("T", "N", j1, j1 - j0 + 1, nt, 1.0_dp, Y, nt, Y(1, j0), nt, &
                   1.0_dp, A(1, j0), nk)
    end do
    !$omp end parallel do
    deallocate(X, Y)

    ! ---- Moebius error s_n(1/2) from the same Gram data (check) ----
    allocate(mu(nmax), spf(nmax), sg(nmax))
    call sieve(nmax, mu, spf)
    sg(1) = E2
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

    ! ---- Cholesky A = U^T U and y = U^{-T} b; d_n^2 for every n ----
    call dpotrf("U", nk, A, nk, info)
    nok = nk
    if (info > 0) then
        nok = info - 1
        print '(a,i0,a,i0)', "Cholesky stops at column ", info, ": d_n computed for n <= ", nok + 1
    end if
    allocate(y_(nok), d2(nok + 1))
    y_ = b(1:nok)
    call dtrsv("U", "T", "N", nok, A, nk, y_, 1)
    d2(1) = E2
    do n1 = 1, nok
        d2(n1 + 1) = d2(n1) - y_(n1)**2
    end do

    call system_clock(clock1)
    ttotal = real(clock1 - clock0, dp) / real(rate, dp)

    ! ---- output: "n d_n^2 s_n" separated by one space ----
    open(newunit=u, file=trim(out_file), status="replace", action="write")
    write(u, '(a)') "# d_n(1/2)^2 (optimal error) and s_n(1/2) (Moebius error) from the Gram matrix"
    write(u, '(a,i0,a,es9.2,a,es9.2,a,i0,a,i0,a,f10.2,a)') "# nodes = ", nt, ", T_fine = ", &
        T_fine, ", T_max = ", T_max, ", nmax = ", nmax, ", stable up to n = ", nok + 1, &
        ", wall time = ", ttotal, " s"
    write(u, '(a,es23.15)') "# ||E||^2 on the grid = ", E2
    do k = 2, nmax, max(1, (nmax - 2) / 6)
        write(u, '(a,i0,a,es23.15,a,es23.15)') "# <E,G_k> for k = ", k, ": ", b(k - 1), &
            "   exact -log(k)/k = ", -log(real(k, dp)) / real(k, dp)
    end do
    write(u, '(a)') "# n d_n^2 s_n"
    do n = 2, nok + 1
        line = inum(n) // " " // num(d2(n), '(es16.8)') // " " // num(sg(n), '(es16.8)')
        write(u, '(a)') line
    end do
    close(u)

    print '(a,i0,a,f8.6,a,f10.2,a,a)', "done: n <= ", nok + 1, ", ||E||^2 = ", E2, &
        ", wall time ", ttotal, " s  ->  ", trim(out_file)
end program riemann_dn
