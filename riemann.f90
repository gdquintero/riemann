program riemann
    use stdlib_kinds,       only: sp, dp, int64

    implicit none

    integer, parameter :: I128 = selected_int_kind(38)
    real(dp), parameter :: pi = 3.1415926535897932
    real(dp) :: alpha_saved
    integer(selected_int_kind(38)) :: n_saved

    contains

    subroutine norm_hardy(alpha,n,norm,ier_hardy)
        use stdlib_kinds, only: dp
        use quadpack,     only: dqagi
        use quadpack,     only: dqag

        implicit none

        integer, parameter :: I128 = selected_int_kind(38)

        real(dp),      intent(in)  :: alpha
        integer(I128), intent(in)  :: n
        real(dp),      intent(out) :: norm
        integer,       intent(out) :: ier_hardy

        real(dp), parameter :: bound  = 0.0_dp
        integer,  parameter :: inf    = 1
        real(dp), parameter :: epsabs = 1.0e-10_dp
        real(dp), parameter :: epsrel = 1.0e-10_dp

        integer,  parameter :: limit  = 200
        integer,  parameter :: lenw   = 4*limit

        real(dp) :: abserr, work(lenw)
        integer  :: neval, last, iwork(limit)

        real(dp), parameter :: a = 0.0_dp, b = 100.0_dp
        integer,  parameter :: key = 6

        ! Guardar parámetros para el integrando
        alpha_saved = alpha
        n_saved     = n

        call dqagi(hardy_f, bound, inf, epsabs, epsrel, norm, abserr, neval, ier_hardy, &
                    limit, lenw, last, iwork, work)

    end subroutine norm_hardy


    real(dp) function hardy_f(t)
        use stdlib_kinds, only: dp

        implicit none

        real(dp), intent(in) :: t
        complex(dp) :: Gn

        call G(alpha_saved, t, n_saved, Gn)     ! tu G devuelve complejo
        hardy_f = real(Gn*conjg(Gn), dp) / pi   ! pi debe estar en el host del programa
    end function hardy_f

    subroutine G(alpha, t, n, Gn)
        use stdlib_kinds, only: dp

        implicit none

        integer, parameter :: I128 = selected_int_kind(38)
        real(dp),      intent(in)  :: alpha, t
        integer(I128), intent(in)  :: n
        complex(dp),   intent(out) :: Gn

        complex(dp) :: s, zeta, Hk
        integer(I128) :: k
        integer :: mu

        s = cmplx(alpha, t, kind=dp)

        ! zeta(s) SOLO depende de s: calcular una vez
        call zeta_function(s, zeta)

        Gn = (0.0_dp, 0.0_dp)

        do k = 2, n
            call mobius(k, mu)
            if (mu == 0) cycle

            call Hk_eval(s, k, zeta, Hk)
            Gn = Gn + (real(mu, dp) / real(k, dp)) * Hk
        end do

        Gn = Gn + 1.0_dp / s
    end subroutine G

    subroutine Hk_eval(s, k, zeta, Hk)
        use stdlib_kinds, only: dp

        implicit none

        integer, parameter :: I128 = selected_int_kind(38)
        complex(dp),   intent(in)  :: s
        integer(I128), intent(in)  :: k
        complex(dp),   intent(in)  :: zeta
        complex(dp),   intent(out) :: Hk

        real(dp)    :: logk
        complex(dp) :: kpow, factor

        ! k^(1-s) = exp((1-s)*log(k))
        logk = log(real(k, dp))
        kpow = exp( (1.0_dp - s) * logk )

        factor = (1.0_dp - kpow)

        Hk = factor * (zeta / s)
    end subroutine Hk_eval

    subroutine zeta_floor(s,n,zeta)
        use stdlib_kinds, only: dp

        implicit none

        complex(dp), intent(in)  :: s
        integer,     intent(in)  :: n
        complex(dp), intent(out) :: zeta

        integer :: k
        real(dp) :: logk, logkp1
        complex(dp) :: ks, kp1s, k1ms, kp11ms
        complex(dp) :: term, one_minus_s

        if (abs(s - cmplx(1.0_dp, 0.0_dp, kind=dp)) < 1.0e-14_dp) then
            zeta = cmplx(0.0_dp, 0.0_dp, kind=dp)
            return
        end if

        if (abs(s) < 1.0e-14_dp .or. abs(1.0_dp - s) < 1.0e-14_dp) then
            zeta = cmplx(0.0_dp, 0.0_dp, kind=dp)
            return
        end if

        ! Parte racional: s/(s-1)
        zeta = s/(s - 1.0_dp)

        ! Suma truncada en k=1..n
        do k = 1, n
            logk   = log(real(k,   dp))
            logkp1 = log(real(k+1, dp))

            ! k^{-s}, (k+1)^{-s}
            ks   = exp((-s) * logk)
            kp1s = exp((-s) * logkp1)

            ! k^{1-s}, (k+1)^{1-s}
            k1ms   = exp((1.0_dp - s) * logk)
            kp11ms = exp((1.0_dp - s) * logkp1)

            term = real(k,dp) * (ks - kp1s) / s  -  (kp11ms - k1ms) / (1.0_dp - s)
            zeta = zeta + s * term
        end do

    end subroutine zeta_floor

    subroutine zeta_function(s,zeta)
        use stdlib_kinds, only: dp

        implicit none

        complex(dp),    intent(in)  :: s
        complex(dp),    intent(out) :: zeta

        integer :: n,n_max,n_used
        real(dp) :: tol,logk,logkp1
        complex(dp) :: ks,kp1s,k1ms,kp11ms,term,zeta1,zeta2

        n = 1
        n_max = 1000000
        tol = 1.0e-10_dp

        call zeta_floor(s,n,zeta1)

        do
            if (2 * n > n_max) then
                zeta = zeta1
                n_used = n
                return
            endif

            call zeta_floor(s,2*n,zeta2)

            if (abs(zeta2 - zeta1) < tol) then
                zeta = zeta2
                n_used = n
                return
            endif

            n = 2 * n
            zeta1 = zeta2
        enddo
    end subroutine zeta_function

    subroutine mobius(n, mu)
        implicit none

        integer, parameter :: I128 = selected_int_kind(38)
        integer(I128), intent(in)  :: n
        integer,       intent(out) :: mu

        integer(I128) :: m, p
        integer :: k

        if (I128 < 0) then
            mu = 0
            return
        end if

        if (n <= 0_I128) then
            mu = 0
            return
        end if
        if (n == 1_I128) then
            mu = 1
            return
        end if

        m = n
        k = 0

        ! Factor 2
        if (mod(m, 2_I128) == 0_I128) then
            m = m / 2_I128
            if (mod(m, 2_I128) == 0_I128) then
            mu = 0
            return
            end if
            k = k + 1
        end if

        ! Factor 3
        if (mod(m, 3_I128) == 0_I128) then
            m = m / 3_I128
            if (mod(m, 3_I128) == 0_I128) then
            mu = 0
            return
            end if
            k = k + 1
        end if

        ! División de prueba 6t ± 1
        ! Condición segura: p <= m/p (evita overflow en p*p)
        p = 5_I128
        do while (p <= m / p)

            if (mod(m, p) == 0_I128) then
            m = m / p
            if (mod(m, p) == 0_I128) then
                mu = 0
                return
            end if
            k = k + 1
            end if

            if (mod(m, p + 2_I128) == 0_I128) then
            m = m / (p + 2_I128)
            if (mod(m, p + 2_I128) == 0_I128) then
                mu = 0
                return
            end if
            k = k + 1
            end if

            p = p + 6_I128
        end do

        ! Si queda un primo > 1, cuenta como factor distinto
        if (m > 1_I128) k = k + 1

        if (mod(k, 2) == 0) then
            mu = 1
        else
            mu = -1
        end if
    end subroutine mobius
    
end program riemann
