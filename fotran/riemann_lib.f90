! Common routines of riemann_tol and riemann_traj.
!
! Both programs evaluate s_n(alpha) = || S_n - E ||^2_alpha, S_n = sum_{k=2}^n mu(k) G_k,
! through the exact identity (inverse Mellin transform, u = 1/x)
!
!     s_n(alpha) = sum_{m >= 1} w_m(alpha) R_n(m)^2 ,
!     R_n(m)     = 1 + sum_{k=2}^{n} mu(k) (m mod k) / k ,
!     w_m(alpha) = (m^(-2 alpha) - (m+1)^(-2 alpha)) / (2 alpha) .
!
! The sum over m is split into m <= m0 (exact weights), geometric blocks m0 < m <= mmax
! (w expanded to second order around the block centre, three moments of R^2 per block)
! and the tail m > mmax (R^2 replaced by its exact mean over a period, see tail_mean).

module riemann_lib
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: iso_c_binding,   only: c_double
    implicit none
    private
    public :: dp, pow_diff, sieve, squarefree_divisors, geometric_blocks, weights, &
              tail_mean, num, inum

    interface
        pure function c_expm1(x) bind(C, name="expm1") result(y)
            import :: c_double
            real(c_double), value :: x
            real(c_double)        :: y
        end function c_expm1

        pure function c_log1p(x) bind(C, name="log1p") result(y)
            import :: c_double
            real(c_double), value :: x
            real(c_double)        :: y
        end function c_log1p
    end interface

contains

    ! x^(-p) - (x+1)^(-p), without cancellation for large x
    elemental real(dp) function pow_diff(x, p)
        real(dp), intent(in) :: x, p

        pow_diff = -x**(-p) * c_expm1(-p * c_log1p(1.0_dp / x))
    end function pow_diff

    ! Linear sieve: Moebius function and smallest prime factor up to nmax
    subroutine sieve(nmax, mu, spf)
        integer, intent(in)  :: nmax
        integer, intent(out) :: mu(nmax), spf(nmax)

        integer, allocatable :: primes(:)
        integer :: i, j, p, np

        allocate(primes(nmax))
        np  = 0
        spf = 0
        mu(1) = 1
        do i = 2, nmax
            if (spf(i) == 0) then
                spf(i) = i
                mu(i)  = -1
                np = np + 1
                primes(np) = i
            end if
            do j = 1, np
                p = primes(j)
                if (p > spf(i) .or. p > nmax / i) exit
                spf(i*p) = p
                if (p == spf(i)) then
                    mu(i*p) = 0
                else
                    mu(i*p) = -mu(i)
                end if
            end do
        end do
    end subroutine sieve

    ! Divisors d >= 2 of a squarefree n, with J_2(d) = prod_{p|d} (p^2 - 1)
    subroutine squarefree_divisors(n, spf, nd, d, j2)
        integer,  intent(in)  :: n, spf(:)
        integer,  intent(out) :: nd, d(:)
        real(dp), intent(out) :: j2(:)

        integer  :: p(16), np, q, mask, i, dd
        real(dp) :: jj

        np = 0
        q  = n
        do while (q > 1)
            np = np + 1
            p(np) = spf(q)
            q = q / p(np)
        end do

        nd = 0
        do mask = 1, 2**np - 1
            dd = 1
            jj = 1.0_dp
            do i = 1, np
                if (btest(mask, i - 1)) then
                    dd = dd * p(i)
                    jj = jj * (real(p(i), dp)**2 - 1.0_dp)
                end if
            end do
            nd = nd + 1
            d(nd)  = dd
            j2(nd) = jj
        end do
    end subroutine squarefree_divisors

    ! Blocks [bstart(i), bstart(i+1)) covering m0+1 .. mmax, of relative width h and at
    ! most cap terms each, with their centres
    subroutine geometric_blocks(m0, mmax, h, cap, nb, bstart, centre)
        integer,  intent(in)  :: m0, mmax, cap
        real(dp), intent(in)  :: h
        integer,  intent(out) :: nb
        integer,  allocatable, intent(out) :: bstart(:)
        real(dp), allocatable, intent(out) :: centre(:)

        integer :: b, ib

        nb = 0
        b  = m0 + 1
        do while (b <= mmax)
            nb = nb + 1
            b  = next_start(b)
        end do
        allocate(bstart(nb + 1), centre(nb))
        bstart(1) = m0 + 1
        do ib = 1, nb
            bstart(ib + 1) = next_start(bstart(ib))
            centre(ib) = 0.5_dp * real(bstart(ib) + bstart(ib + 1) - 1, dp)
        end do

    contains

        integer function next_start(b)
            integer, intent(in) :: b

            next_start = min(mmax + 1, b + cap, max(b + 1, int(real(b, dp) * (1.0_dp + h))))
        end function next_start
    end subroutine geometric_blocks

    ! Exact weights for m <= m0, Taylor coefficients (w, w', w''/2) at the block centres
    ! and the total weight of the tail m > mmax, for every alpha
    subroutine weights(alpha, m0, centre, mmax, wexact, wb, tailw)
        real(dp), intent(in)  :: alpha(:), centre(:)
        integer,  intent(in)  :: m0, mmax
        real(dp), intent(out) :: wexact(:,:), wb(:,:,:), tailw(:)

        integer  :: ia, m
        real(dp) :: p

        do ia = 1, size(alpha)
            p = 2.0_dp * alpha(ia)
            wexact(:, ia) = pow_diff([(real(m, dp), m = 1, m0)], p) / p
            wb(1, :, ia)  = pow_diff(centre, p) / p                               ! w
            wb(2, :, ia)  = -pow_diff(centre, p + 1.0_dp)                         ! w'
            wb(3, :, ia)  = 0.5_dp * (p + 1.0_dp) * pow_diff(centre, p + 2.0_dp)  ! w''/2
            tailw(ia)     = real(mmax + 1, dp)**(-p) / p
        end do
    end subroutine weights

    ! Add the term k = n (mu(n) /= 0) to the exact period mean of R^2,
    !   E[R^2] = E[R]^2 + (1/12) sum_{d>=2} J_2(d) T_d^2,  T_d = sum_{2<=k<=n, d|k} mu(k)/k,
    ! kept as ER = E[R], Q = sum_d J_2(d) T_d^2 and T(d).  Returns the updated E[R^2].
    real(dp) function tail_mean(n, mu_n, spf, ER, Q, T)
        integer,  intent(in)    :: n, mu_n, spf(:)
        real(dp), intent(inout) :: ER, Q, T(:)

        integer  :: dv(512), nd, i
        real(dp) :: j2(512), c

        c  = real(mu_n, dp) / real(n, dp)
        ER = ER + 0.5_dp * c * real(n - 1, dp)
        call squarefree_divisors(n, spf, nd, dv, j2)
        do i = 1, nd
            Q = Q + j2(i) * (2.0_dp * T(dv(i)) * c + c * c)
            T(dv(i)) = T(dv(i)) + c
        end do
        tail_mean = ER * ER + Q / 12.0_dp
    end function tail_mean

    ! Text of x in format f without surrounding blanks
    function num(x, f) result(str)
        real(dp),     intent(in)  :: x
        character(*), intent(in)  :: f
        character(:), allocatable :: str

        character(len=40) :: buf

        write(buf, f) x
        str = trim(adjustl(buf))
    end function num

    function inum(k) result(str)
        integer, intent(in)       :: k
        character(:), allocatable :: str

        character(len=20) :: buf

        write(buf, '(i0)') k
        str = trim(buf)
    end function inum

end module riemann_lib
