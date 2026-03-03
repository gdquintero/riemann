program riemann
    use stdlib_kinds, only: sp, dp, int64

    implicit none

    integer, parameter :: I128 = selected_int_kind(38)
    integer :: nodes
    integer(I128) :: k
    real(dp), allocatable :: x(:), w(:),interval(:)
    complex(dp) :: s,hk
    real(dp) :: alpha, t

    nodes = 256
    allocate(x(nodes), w(nodes),interval(2))

    alpha = 2.0_dp
    t = 3.0_dp
    k = 100
    
    call H(alpha,t,k,x,w,nodes,interval,hk)

    print*, pi

    contains

    subroutine zeta_function(s,x,w,nodes,interval,zeta)
        use stdlib_kinds, only: dp
        use stdlib_quadrature, only: gauss_legendre
        use stdlib_specialfunctions_gamma, only: gamma
        implicit none

        integer,        intent(in) :: nodes
        complex(dp),    intent(in) :: s
        real(dp),       intent(inout) :: x(:), w(:), interval(:)
        complex(dp),    intent(out) :: zeta
        real(dp), parameter :: pi = 3.1415926535897932

        complex(dp) :: a,xi

        interval = [0.0_dp,1.0_dp]
        a = 0.5_dp * s

        call gauss_legendre(x,w,interval)
        call xi_quad(s,x,w,nodes,xi)

        zeta = exp( a * log(pi)) * xi / gamma(a)

    end subroutine zeta_function

    subroutine xi_quad(s,x,w,nodes,xi)
        use stdlib_kinds, only: dp
        implicit none

        integer,        intent(in) :: nodes
        complex(dp),    intent(in)  :: s
        real(dp),       intent(in)  :: x(:),w(:)
        complex(dp),    intent(out) :: xi

        complex(dp) :: a,b,quad,term1,term2
        real(dp) :: u,jac,psi_u,logu
        integer :: i

        a = -0.5_dp * (s + 1.0_dp)
        b =  0.5_dp * s - 1.0_dp

        quad = (0.0_dp, 0.0_dp)

        do i = 1, nodes
            u = 1.0_dp/(1.0_dp - x(i))
            jac = 1.0_dp/(1.0_dp - x(i))**2  
            logu = log(u)

            term1 = exp(a * logu)    ! u^(-s/2-1/2)
            term2 = exp(b * logu)    ! u^( s/2-1)

            call psi(u,psi_u)

            quad = quad + (term1 + term2) * psi_u * jac * w(i)
        end do

        xi = 1.0_dp/(s - 1.0_dp) - 1.0_dp/s + quad
    end subroutine xi_quad

    subroutine norm_hardy(alpha,x,w,nodes,res)
        use stdlib_kinds, only: dp
        implicit none

        real(dp),    intent(in)  :: alpha
        real(dp),    intent(in)  :: x(:), w(:)
        integer,     intent(in)  :: nodes
        real(dp),    intent(out) :: res

        integer :: i
        real(dp) :: t, jac
        complex(dp) :: s, Fva
        real(dp), parameter :: pi = 3.1415926535897932



    end subroutine norm_hardy

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

    subroutine psi(u,res)
        use stdlib_kinds, only: dp
        implicit none

        real(dp),   intent(in) :: u
        real(dp),   intent(out) :: res
        
        real(dp), parameter :: tol = 1.d-8, pi = 4.d0 * atan(1.d0) 
        real(dp) :: term
        integer :: n, n_max

        res = 0.d0
        n = 1
        n_max = 10000
        term = exp(-pi * n**2 * u)

        do while (term .gt. tol .and. n .le. n_max)
            res = res + term
            n = n + 1
            term = exp(-pi * n**2 * u)
        enddo

    end subroutine psi

    subroutine H(alpha,t,k,x,w,nodes,interval,hk)
        use stdlib_kinds, only: dp
        implicit none

        integer, parameter :: I128 = selected_int_kind(38)
        integer(I128),  intent(in)  :: k
        integer,        intent(in) :: nodes
        real(dp),       intent(in) :: alpha,t
        complex(dp),    intent(out) :: hk
        real(dp),       intent(inout) :: x(:), w(:), interval(:)

        complex(dp) :: s,zeta,a

        s = cmplx(alpha,t,kind=dp)
        a = 1.0_dp - exp((1.0_dp - s) * log(real(k, dp)))

        call zeta_function(s,x,w,nodes,interval,zeta)

        hk = a * zeta / s
        
    end subroutine H

    ! subroutine hardy_norm_Gn()
    
end program riemann
