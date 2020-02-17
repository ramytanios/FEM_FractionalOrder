#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

/// Compile-time parameters 
constexpr double theta = 0.5; 

constexpr double tol = 1e-10; 

constexpr double beta = 0.5; /// = 1-\beta = 1-\alpha/2
constexpr bool testbool = (beta==0); /// to distinguish between fractional
				     /// and non-fractional cases

static_assert(theta>=0 && theta<=1, "theta is out of range!"); 
static_assert(beta>=0 && beta<1, "beta is out of range!"); 

#endif 
