! Tolerance experiment of riemann_0.01.py in Fortran: for each alpha, the smallest n with
! s_n(alpha) = || S_n - E ||^2_alpha <= tol, where S_n = sum_{k=2}^n mu(k) G_k.
!
! Instead of integrating |S_n - E|^2 over Re s = alpha (scipy quad + zeta), the error is
! evaluated exactly through the inverse Mellin transform and the change u = 1/x.  On
! [m, m+1) the function whose Mellin transform is E - S_n equals
!
!     R_n(m) = 1 + sum_{k=2}^{n} mu(k) (m mod k) / k ,
!
! so that, with w_m(alpha) = (m^(-2 alpha) - (m+1)^(-2 alpha)) / (2 alpha),
!
!     s_n(alpha) = sum_{m >= 1} w_m(alpha) R_n(m)^2 .
!
! No zeta function and no quadrature in t are needed.  The sum over m is split into
!   * 1 <= m <= m0     : summed exactly, one weight per m and alpha;
!   * m0 < m <= mmax   : geometric blocks of relative width h; w is expanded to second
!                        order around the block centre, so each alpha only needs three
!                        moments of R^2 per block;
!   * m > mmax         : R^2 is replaced by its exact mean over a period,
!                        E[R^2] = E[R]^2 + (1/12) sum_{d>=2} J_2(d) T_d^2,
!                        T_d = sum_{2<=k<=n, d|k} mu(k)/k.
! R is updated in place when n grows, so each new term costs one pass over m, shared by
! all values of alpha.  The pass is memory bound: more than ~8 threads does not help.
!
! Convention: n is the last index of S_n.  The Python scripts (np.arange(2, n)) report
! n + 1 for the same partial sum.
!
! Usage:  ./riemann_tol [file.nml]      (namelist &params, e.g. riemann_0.01.nml)

program riemann_tol
    use riemann_lib
    use, intrinsic :: iso_fortran_env, only: int64
    !$ use omp_lib, only: omp_set_num_threads
    implicit none

    ! ---- parameters (namelist &params); defaults reproduce riemann_0.01.py ----
    real(dp) :: alpha_min, alpha_max, h, tol(8)
    integer  :: npts, nmax, mmax, m0, traj_every, nthreads
    logical  :: full
    character(len=256) :: out_file, traj_file, nml_file
    namelist /params/ alpha_min, alpha_max, npts, tol, nmax, mmax, m0, h, full, &
                      traj_every, out_file, traj_file, nthreads

    ! ---- work arrays ----
    integer,  allocatable :: mu(:), spf(:), bstart(:), nfirst(:,:), nlast(:,:)
    real(dp), allocatable :: alpha(:), R(:), R2(:), wexact(:,:), wb(:,:,:), mom(:,:)
    real(dp), allocatable :: centre(:), tailw(:), s(:), sfirst(:,:), tfirst(:,:), T(:)
    integer  :: na, ntol, nb, ib, ia, it, i, m, n, rem, u, utraj, n_end
    integer(int64) :: clock0, clock1, rate
    real(dp) :: c, t2, dm, a0, a1, a2, ER, Q, ER2, ttotal
    character(len=32) :: fmt
    character(len=:), allocatable :: line

    alpha_min  = 0.5_dp
    alpha_max  = 0.995_dp
    npts       = 100
    tol        = -1.0_dp
    tol(1)     = 1.0e-2_dp
    nmax       = 10000
    mmax       = 10000000
    m0         = 10000
    h          = 1.0e-3_dp
    full       = .false.
    traj_every = 0
    out_file   = "output_tol_0.01.txt"
    traj_file  = "traj_f90.txt"
    nthreads   = 8

    if (command_argument_count() >= 1) then
        call get_command_argument(1, nml_file)
        open(newunit=u, file=trim(nml_file), status="old", action="read")
        read(u, nml=params)
        close(u)
    end if

    ntol = count(tol > 0.0_dp)
    tol(1:ntol) = pack(tol, tol > 0.0_dp)
    if (ntol == 0) error stop "riemann_tol: no positive tolerance given"
    if (nmax < 2) error stop "riemann_tol: nmax must be >= 2"
    m0 = min(max(m0, 1), mmax)
    if (mmax < 100 * nmax) print '(a)', "warning: mmax < 100*nmax, the tail estimate may be poor"
    !$ if (nthreads > 0) call omp_set_num_threads(nthreads)

    call system_clock(clock0, rate)

    ! ---- alpha grid (np.linspace) ----
    na = npts
    allocate(alpha(na))
    if (na == 1) then
        alpha(1) = alpha_min
    else
        alpha = [(alpha_min + (alpha_max - alpha_min) * real(i - 1, dp) / real(na - 1, dp), i = 1, na)]
    end if

    allocate(mu(nmax), spf(nmax), T(nmax))
    call sieve(nmax, mu, spf)

    ! ---- geometric blocks covering m0+1 .. mmax, and the weights ----
    call geometric_blocks(m0, mmax, h, mmax, nb, bstart, centre)
    allocate(mom(3, nb), wexact(m0, na), wb(3, nb, na), tailw(na))
    call weights(alpha, m0, centre, mmax, wexact, wb, tailw)

    ! ---- state for n = 1 (S_1 = 0, residual R = 1) ----
    allocate(R(mmax), R2(m0), s(na))
    R  = 1.0_dp
    T  = 0.0_dp
    ER = 1.0_dp
    Q  = 0.0_dp
    s  = 1.0_dp / (2.0_dp * alpha)

    allocate(nfirst(ntol, na), nlast(ntol, na), sfirst(ntol, na), tfirst(ntol, na))
    nfirst = 0
    nlast  = 1
    sfirst = 0.0_dp
    tfirst = 0.0_dp

    if (traj_every > 0) then
        open(newunit=utraj, file=trim(traj_file), status="replace", action="write")
        write(utraj, '(a)', advance="no") "#        n"
        do ia = 1, na
            write(utraj, '(a,f8.4)', advance="no") "    alpha=", alpha(ia)
        end do
        write(utraj, *)
        write(fmt, '(a,i0,a)') '(i10,', na, 'es18.10)'
    end if

    n_end = nmax
    do n = 2, nmax
        if (mu(n) /= 0) then
            c = real(mu(n), dp) / real(n, dp)

            ! R(m) += mu(n) (m mod n)/n, and moments of R^2 per block
            !$omp parallel do
            do m = 1, m0
                R(m) = R(m) + c * real(mod(m, n), dp)
            end do
            !$omp end parallel do

            !$omp parallel do schedule(dynamic, 4) private(m, rem, t2, dm, a0, a1, a2)
            do ib = 1, nb
                rem = mod(bstart(ib), n)
                a0  = 0.0_dp
                a1  = 0.0_dp
                a2  = 0.0_dp
                do m = bstart(ib), bstart(ib + 1) - 1
                    R(m) = R(m) + c * real(rem, dp)
                    t2 = R(m) * R(m)
                    dm = real(m, dp) - centre(ib)
                    a0 = a0 + t2
                    a1 = a1 + dm * t2
                    a2 = a2 + dm * dm * t2
                    rem = rem + 1
                    if (rem == n) rem = 0
                end do
                mom(1, ib) = a0
                mom(2, ib) = a1
                mom(3, ib) = a2
            end do
            !$omp end parallel do

            ! exact period mean of R^2 for the tail m > mmax
            ER2 = tail_mean(n, mu(n), spf, ER, Q, T)

            R2 = R(1:m0)**2
            !$omp parallel do
            do ia = 1, na
                s(ia) = dot_product(wexact(:, ia), R2) + sum(wb(:, :, ia) * mom) &
                        + tailw(ia) * ER2
            end do
            !$omp end parallel do
        end if

        do ia = 1, na
            do it = 1, ntol
                if (s(ia) <= tol(it)) then
                    if (nfirst(it, ia) == 0) then
                        call system_clock(clock1)
                        nfirst(it, ia) = n
                        sfirst(it, ia) = s(ia)
                        tfirst(it, ia) = real(clock1 - clock0, dp) / real(rate, dp)
                    end if
                else
                    nlast(it, ia) = n
                end if
            end do
        end do

        if (traj_every > 0) then
            if (mod(n, traj_every) == 0) write(utraj, fmt) n, s
        end if

        if (mod(n, max(1, nmax / 10)) == 0) then
            call system_clock(clock1)
            print '(a,i0,a,f9.1,a)', "n = ", n, "   (", real(clock1 - clock0, dp) / real(rate, dp), " s)"
        end if

        if (.not. full .and. all(nfirst > 0)) then
            n_end = n
            exit
        end if
    end do

    if (traj_every > 0) close(utraj)
    call system_clock(clock1)
    ttotal = real(clock1 - clock0, dp) / real(rate, dp)

    ! ---- output: one line per alpha, values separated by one space.  With a single tol
    !      the first four columns are those of the Python scripts: alpha s n t ----
    open(newunit=u, file=trim(out_file), status="replace", action="write")
    write(u, '(a)') "# s_n(alpha) = ||S_n - E||^2_alpha, S_n = sum_{k=2}^n mu(k) G_k, exact (no zeta)"
    write(u, '(a,i0,a,i0,a,i0,a,es9.2,a,l1,a,f10.2,a)') "# nmax = ", nmax, ", mmax = ", mmax, &
        ", m0 = ", m0, ", h = ", h, ", full = ", full, ", wall time = ", ttotal, " s"
    write(u, '(a)') "# per tol: s n t reached n_last, where n = smallest n with s_n <= tol (n_end if"
    write(u, '(a)') "# not reached), s = s_n at that n, t = elapsed seconds, reached = 1/0 and"
    write(u, '(a)') "# n_last = largest n <= n_end with s_n > tol.  n is the last index of S_n;"
    write(u, '(a)') "# the Python scripts reported n + 1 for the same partial sum."
    line = "# alpha"
    do it = 1, ntol
        line = line // " [tol=" // num(tol(it), '(es8.1)') // "] s n t reached n_last"
    end do
    write(u, '(a)') line // " n_end s_end"
    do ia = 1, na
        line = num(alpha(ia), '(f8.4)')
        do it = 1, ntol
            if (nfirst(it, ia) > 0) then
                line = line // " " // num(sfirst(it, ia), '(es16.8)') // " " // inum(nfirst(it, ia)) &
                       // " " // num(tfirst(it, ia), '(f12.3)') // " 1 " // inum(nlast(it, ia))
            else
                line = line // " " // num(s(ia), '(es16.8)') // " " // inum(n_end) &
                       // " " // num(ttotal, '(f12.3)') // " 0 " // inum(nlast(it, ia))
            end if
        end do
        write(u, '(a)') line // " " // inum(n_end) // " " // num(s(ia), '(es16.8)')
    end do
    close(u)

    print '(a,f10.2,a,a)', "done in ", ttotal, " s  ->  ", trim(out_file)

end program riemann_tol
