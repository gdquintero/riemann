program riemann
    use stdlib_kinds,                   only: sp, dp, int64
    use stdlib_quadrature,              only: gauss_legendre

    implicit none

    integer, parameter :: I128 = selected_int_kind(38)
    real(dp), parameter :: pi = 3.1415926535897932
    integer :: nodes
    integer(I128) :: n
    real(dp), allocatable :: x(:), w(:),interval(:)
    real(dp) :: alpha, norm

    nodes = 10
    allocate(x(nodes), w(nodes),interval(2))

    interval = [0.0_dp,1.0_dp]
    call gauss_legendre(x,w,interval)

    n = 2
    alpha = 0.8_dp

    call norm_hardy(alpha,n,x,w,nodes,norm)

    ! print*, norm
    

    contains

    subroutine norm_hardy(alpha,n,norm,ier)
        use stdlib_kinds, only: dp
        use quadpack,     only: dqagi

        implicit none

        integer, parameter :: I128 = selected_int_kind(38)

        real(dp),       intent(in)  :: alpha
        integer(I128),  intent(in)  :: n
        real(dp),       intent(out) :: norm
        integer,        intent(out) :: ier

        real(dp),       parameter   :: bound = 0.0_dp, epsabs = 1.0e-10_dp, epsrel = 1.0e-10_dp
        integer,        parameter   :: inf   = 1, limit = 200
        integer,        parameter   :: lenw  = 4*limit

        real(dp)                    :: abserr, work(lenw)
        integer                     :: neval, last, iwork(limit)

        call dqagi(f, bound, inf, epsabs, epsrel, norm, abserr, neval, ier, &
             limit, lenw, last, iwork, work)

        contains

        real(dp) function f(t)
            real(dp), intent(in) :: t
            complex(dp) :: Gn
            call G(alpha, t, n, Gn)   ! tu G: devuelve Gn complejo
            f = real(Gn*conjg(Gn), dp) / pi
        end function f

        ! integer :: i
        ! real(dp) :: t, jac
        ! complex(dp) :: s, Gn

        ! res = 0.0_dp

        ! do i = 1, nodes
        !     t   = x(i) / (1.0_dp - x(i))          ! (0,1)->(0,inf)
        !     jac = 1.0_dp / (1.0_dp - x(i))**2     ! dt/dx

        !     s = cmplx(alpha, t, kind=dp)

        !     call G(alpha,t,n,x,w,nodes,Gn)                      
        !     res = res + real(Gn*conjg(Gn), dp) * jac * w(i)
        !     ! En norm_hardy, dentro del loop:

        !     ! print*, "i=", i, "t=", t
        !     ! call G(alpha, t, n, x, w, nodes, Gn)
        !     ! print*, "Gn=", Gn
        ! end do

        ! res = res / pi

    end subroutine norm_hardy

    subroutine G(alpha,t,n,res)
        use stdlib_kinds, only: dp
        implicit none

        integer, parameter :: I128 = selected_int_kind(38)
        real(dp),       intent(in) :: alpha,t
        integer(I128),  intent(in) :: n
        complex(dp),    intent(out):: res

        integer(I128) :: k
        integer :: muk
        complex(dp) :: Hk

        res = (0.0_dp,0.0_dp)

        do k = 2, n
            call mobius(k,muk)
            call H(alpha,t,k,Hk)
            res = res + (real(muk,dp) / real(k,dp) ) * Hk
        enddo

        res = res + 1.0_dp / cmplx(alpha,t,kind=dp)

    end subroutine G

    subroutine H(alpha,t,k,Hk)
        use stdlib_kinds, only: dp
        implicit none

        integer, parameter :: I128 = selected_int_kind(38)
        integer(I128),  intent(in)  :: k
        real(dp),       intent(in) :: alpha,t
        complex(dp),    intent(out) :: Hk

        complex(dp) :: s,zeta,a

        s = cmplx(alpha,t,kind=dp)
        a = 1.0_dp - exp((1.0_dp - s) * log(real(k, dp)))

        call zeta_function(s,zeta)

        hk = a * zeta / s
        
    end subroutine H

    subroutine zeta_function(s,zeta)
        use stdlib_kinds,                   only: dp
        use stdlib_specialfunctions_gamma,  only: gamma

        implicit none

        complex(dp),    intent(in)  :: s
        complex(dp),    intent(out) :: zeta

        complex(dp) :: a,xi

        a = 0.5_dp * s
        
        call xi_quad(s,xi)

        zeta = exp( a * log(pi)) * xi / gamma(a)

    end subroutine zeta_function

    subroutine xi_quad(s,xi)
        use stdlib_kinds, only: dp

        implicit none

        complex(dp),    intent(in)  :: s
        complex(dp),    intent(out) :: xi

        complex(dp) :: a,b,quad,term1,term2
        real(dp)    :: u,jac,psi_u,logu
        integer     :: i

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

    subroutine psi(u,res)
        use stdlib_kinds, only: dp
        implicit none

        real(dp),   intent(in) :: u
        real(dp),   intent(out) :: res
        
        real(dp), parameter :: tol = 1.d-8
        real(dp) :: term
        integer :: n, n_max

        res = 0.d0
        n = 1
        n_max = 10000
        term = exp(-pi * n**2 * u)

        do while (term .gt. tol .and. n .le. n_max)
            res = res + term
            n = n + 1
            term = exp(-pi * n*n * u)
        enddo

    end subroutine psi

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
