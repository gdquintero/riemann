! Experiments of riemann_1.py ... riemann_20.py (launched together by run_all.py) in
! Fortran: s_n(alpha) = || S_n - E ||^2_alpha at n = n_start, n_start + n_step, ..., nmax
! (or, with per_decade > 0, at per_decade log-spaced values of n per decade) for a list
! of alphas, all in one run.
!
! The error is evaluated exactly, without zeta or quadrature (see riemann_lib.f90):
! s_n(alpha) = sum_m w_m(alpha) R_n(m)^2 with R_n(m) = 1 + sum_{k=2}^n mu(k) (m mod k)/k.
! Instead of sweeping all m for every n (memory bound), the range of m is cut into
! tiles that fit in cache and every tile runs through all n on its own: R lives in the
! tile, the work is compute bound and scales with the number of cores.
!
! Convention: n is the last index of S_n.  The Python scripts (np.arange(2, n)) used
! S_{n-1}; for n = 1000 j both coincide because mu(1000 j) = 0.
!
! Usage:  ./riemann_traj [file.nml]     (namelist &params, e.g. riemann_traj.nml)

program riemann_traj
    use riemann_lib
    use, intrinsic :: iso_fortran_env, only: int64
    !$ use omp_lib, only: omp_set_num_threads
    implicit none

    integer, parameter :: exact_tile = 4096

    ! ---- parameters (namelist &params); defaults reproduce the runs of the paper ----
    real(dp) :: alphas(64), h
    integer  :: n_start, n_step, per_decade, nmax, mmax, m0, cap, nthreads
    character(len=256) :: out_dir, nml_file
    namelist /params/ alphas, n_start, n_step, per_decade, nmax, mmax, m0, h, cap, &
                      nthreads, out_dir

    ! ---- work arrays ----
    integer,  allocatable :: mu(:), spf(:), smap(:), samples(:), bstart(:), estart(:)
    real(dp), allocatable :: alpha(:), T(:), er2s(:), centre(:), wexact(:,:), wb(:,:,:)
    real(dp), allocatable :: tailw(:), mom(:,:,:), pe(:,:,:), s(:,:)
    integer  :: na, ns, ntop, nb, ne, n, i, j, k, ia, is, ib, ie, u
    integer(int64) :: clock0, clock1, rate
    real(dp) :: ER, Q, ER2, ttotal
    character(len=512) :: fname
    character(len=:), allocatable :: line

    alphas     = -1.0_dp
    alphas(1:20) = [(0.50_dp + 0.01_dp * real(i, dp), i = 0, 9), &
                    (0.90_dp + 0.01_dp * real(i, dp), i = 0, 9)]
    n_start    = 1000
    n_step     = 1000
    per_decade = 0
    nmax       = 100000
    mmax       = 100000000
    m0         = 10000
    h          = 1.0e-3_dp
    cap        = 32768
    nthreads   = 20
    out_dir    = "traj"

    if (command_argument_count() >= 1) then
        call get_command_argument(1, nml_file)
        open(newunit=u, file=trim(nml_file), status="old", action="read")
        read(u, nml=params)
        close(u)
    end if

    na = count(alphas > 0.0_dp)
    alpha = pack(alphas, alphas > 0.0_dp)
    if (na == 0) error stop "riemann_traj: no positive alpha given"
    if (n_start < 2 .or. n_step < 1 .or. nmax < n_start) error stop "riemann_traj: bad n range"
    m0 = min(max(m0, 1), mmax)
    if (mmax < 100 * nmax) print '(a)', "warning: mmax < 100*nmax, the tail estimate may be poor"
    !$ if (nthreads > 0) call omp_set_num_threads(nthreads)

    call system_clock(clock0, rate)

    ! ---- sample points n (arithmetic or log spaced) and Moebius function ----
    if (per_decade > 0) then
        ns = int(real(per_decade, dp) * log10(real(nmax, dp) / real(n_start, dp))) + 1
        allocate(samples(ns + 1))
        k = 0
        do is = 0, ns - 1
            n = nint(real(n_start, dp) * 10.0_dp**(real(is, dp) / real(per_decade, dp)))
            if (k == 0) then
                k = 1
                samples(k) = n
            else if (n > samples(k)) then
                k = k + 1
                samples(k) = n
            end if
        end do
        if (samples(k) < nmax) then
            k = k + 1
            samples(k) = nmax
        end if
        samples = samples(1:k)
        ns = k
    else
        ns = (nmax - n_start) / n_step + 1
        samples = [(n_start + (is - 1) * n_step, is = 1, ns)]
    end if
    ntop = samples(ns)
    allocate(mu(ntop), spf(ntop), smap(ntop), T(ntop), er2s(ns))
    call sieve(ntop, mu, spf)
    smap = 0
    smap(samples) = [(is, is = 1, ns)]

    ! ---- exact period mean of R^2 (tail m > mmax) at the samples ----
    ER  = 1.0_dp
    Q   = 0.0_dp
    ER2 = 1.0_dp
    T   = 0.0_dp
    do n = 2, ntop
        if (mu(n) /= 0) ER2 = tail_mean(n, mu(n), spf, ER, Q, T)
        if (smap(n) > 0) er2s(smap(n)) = ER2
    end do

    ! ---- tiles: exact ones on 1..m0, geometric blocks on m0+1..mmax; weights ----
    ne = (m0 + exact_tile - 1) / exact_tile
    estart = [(1 + (ie - 1) * exact_tile, ie = 1, ne), m0 + 1]
    call geometric_blocks(m0, mmax, h, cap, nb, bstart, centre)
    allocate(wexact(m0, na), wb(3, nb, na), tailw(na))
    call weights(alpha, m0, centre, mmax, wexact, wb, tailw)

    ! ---- every tile runs through all n; largest tiles first ----
    allocate(mom(3, ns, nb), pe(na, ns, ne))
    !$omp parallel do schedule(dynamic, 1) private(ib, ie)
    do j = 1, nb + ne
        if (j <= nb) then
            ib = nb - j + 1
            call sweep(bstart(ib), bstart(ib + 1), centre(ib), .false., mom(:, :, ib))
        else
            ie = j - nb
            call sweep(estart(ie), estart(ie + 1), 0.0_dp, .true., pe(:, :, ie))
        end if
    end do
    !$omp end parallel do

    ! ---- s_n(alpha) = exact part + Taylor blocks + tail ----
    allocate(s(ns, na))
    do ia = 1, na
        do is = 1, ns
            s(is, ia) = sum(pe(ia, is, :)) + sum(wb(:, :, ia) * mom(:, is, :)) &
                        + tailw(ia) * er2s(is)
        end do
    end do

    call system_clock(clock1)
    ttotal = real(clock1 - clock0, dp) / real(rate, dp)

    ! ---- output: one file per alpha, "n alpha s" separated by one space ----
    call execute_command_line("mkdir -p " // trim(out_dir))
    do ia = 1, na
        if (abs(100.0_dp * alpha(ia) - anint(100.0_dp * alpha(ia))) < 1.0e-9_dp) then
            write(fname, '(a,a,i3.3,a)') trim(out_dir), "/output_", nint(100.0_dp * alpha(ia)), ".txt"
        else
            write(fname, '(a,a,i4.4,a)') trim(out_dir), "/output_", nint(1000.0_dp * alpha(ia)), ".txt"
        end if
        open(newunit=u, file=trim(fname), status="replace", action="write")
        write(u, '(a)') "# s_n(alpha) = ||S_n - E||^2_alpha, S_n = sum_{k=2}^n mu(k) G_k, exact (no zeta)"
        write(u, '(a,i0,a,i0,a,es9.2,a,i0,a,f10.2,a)') "# mmax = ", mmax, ", m0 = ", m0, &
            ", h = ", h, ", cap = ", cap, ", wall time (all alphas) = ", ttotal, " s"
        write(u, '(a)') "# n alpha s_n"
        do is = 1, ns
            line = inum(samples(is)) // " " // num(alpha(ia), '(f8.4)') // " " // &
                   num(s(is, ia), '(es16.8)')
            write(u, '(a)') line
        end do
        close(u)
    end do

    print '(a,i0,a,i0,a,f10.2,a,a,a)', "done: ", na, " alphas x ", ns, " values of n in ", &
        ttotal, " s  ->  ", trim(out_dir), "/"

contains

    ! Run R(m) = R_n(m), a <= m < b, from n = 1 to ntop.  At each sample store either the
    ! weighted sums per alpha (exact = .true.) or the three moments of R^2 around c.
    subroutine sweep(a, b, c, exact, res)
        integer,  intent(in)  :: a, b
        real(dp), intent(in)  :: c
        logical,  intent(in)  :: exact
        real(dp), intent(out) :: res(:,:)

        real(dp), allocatable :: Rt(:), d(:), t2(:)
        integer  :: nlen, n, i, i0, i1, base, is, ia
        real(dp) :: cn, off

        nlen = b - a
        allocate(Rt(nlen), t2(nlen))
        Rt = 1.0_dp
        if (.not. exact) d = [(real(a + i - 1, dp) - c, i = 1, nlen)]

        do n = 2, ntop
            if (mu(n) /= 0) then
                cn = real(mu(n), dp) / real(n, dp)
                ! on each run of m with the same floor(m/n), m mod n = m - base
                base = (a / n) * n
                i0 = 1
                do while (i0 <= nlen)
                    i1  = min(nlen, base + n - a)
                    off = real(a - 1 - base, dp)
                    do i = i0, i1
                        Rt(i) = Rt(i) + cn * (real(i, dp) + off)
                    end do
                    i0   = i1 + 1
                    base = base + n
                end do
            end if

            is = smap(n)
            if (is > 0) then
                t2 = Rt * Rt
                if (exact) then
                    do ia = 1, na
                        res(ia, is) = dot_product(wexact(a:b-1, ia), t2)
                    end do
                else
                    res(1, is) = sum(t2)
                    res(2, is) = dot_product(d, t2)
                    res(3, is) = dot_product(d * d, t2)
                end if
            end if
        end do
    end subroutine sweep

end program riemann_traj
