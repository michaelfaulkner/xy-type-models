function get_monte_carlo_error(mean, mean_squared)
use variables
implicit none
double precision :: get_monte_carlo_error, mean, mean_squared
get_monte_carlo_error = ((mean_squared - mean ** 2) / dfloat(no_of_samples - 1)) ** 0.5
end function get_monte_carlo_error
