module modhelpers

    implicit none

    double precision, parameter :: PI = acos( - 1.0d0)

contains

    subroutine helpers_randomspin(s_, radius_)

        double precision, dimension(3), intent(inout) :: s_
        double precision, optional, intent(in) :: radius_

        double precision :: r, theta, phi

        if (present(radius_)) then
            r = radius_
        else
            r = 1.0
        end if

        call random_number(theta)
        call random_number(phi)

        theta = PI * theta
        phi = 2.0d0 * PI * phi

        s_(1) = r * sin(theta) * cos(phi)
        s_(2) = r * sin(theta) * sin(phi)
        s_(3) = r * cos(theta)

    end subroutine helpers_randomspin

    function helpers_distance(a_, b_, dims_)

        double precision, dimension(:), intent(in) :: a_
        double precision, dimension(:), intent(in) :: b_
        double precision, optional, dimension(:), intent(in) &
            :: dims_
        double precision :: helpers_distance

        if (present(dims_)) then
            helpers_distance = sqrt(sum(min( &
                abs(b_ - a_),                &
                dims_ - abs(b_ - a_)) ** 2))
        else
            helpers_distance = sqrt(sum((b_ - a_) ** 2))
        end if

    end function helpers_distance

end module modhelpers
