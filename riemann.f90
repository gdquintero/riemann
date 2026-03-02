program riemann
    use stdlib_kinds, only: sp, dp, int64
    use stdlib_specialfunctions_gamma, only: gamma
    use stdlib_quadrature, only: gauss_legendre
    implicit none

    complex(sp) :: z
    real(kind=8) :: u,res

    u = 1.d0



    contains

    ! subroutine zeta(z,res)
    !     implicit none


    ! end subroutine zeta

    subroutine psi(u,res)
        implicit none

        real(kind=8),   intent(in) :: u
        real(kind=8),   intent(out) :: res
        
        real(kind=8), parameter :: tol = 1.d-8, pi = 4.d0 * atan(1.d0) 
        real(kind=8) :: term
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

    

    ! subroutine H(alpha,t,n)
    !     implicit none



        
    ! end subroutine H
    
end program riemann
