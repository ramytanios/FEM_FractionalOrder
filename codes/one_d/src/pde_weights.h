#ifndef PDE_WEIGHTS_HPP
#define PDE_WEIGHTS_HPP

#include "params.h"

namespace PdeWeights{

  /**
   * @function alpha_x
   *
   * @brief Evaluate the coefficient functor alpha(x), 
            scalar-valued
   *
   * @param x the coordinates of the point to evaluate the function at
   * @return a scalar the value of alpha at 'x'
   *
   */ 
  std::function<double(double)> alpha_x =
    [] (double x)->double{
      (void)x; 
      return 1; 
};
  
  /**
   * @function gamma_x
   *
   * @brief Evaluate the coefficient functor gamma(x), 
            scalar-valued
   *
   * @param x the coordinates of the point to evaluate the function at
   * @return a scalar which is the value of gamma at 'x'
   *
   */ 
  std::function<double(double)> gamma_x = 
    [] (double x)-> double{
      (void)x; 
      return r;
    };
  
  /**
   * @function beta_x
   *
   * @brief Evaluate the coefficient functor beta(x), 
            scalar-valued
   *
   * @param x the coordinates of the point to evaluate the function at
   * @return a scalar the value of beta at 'x'
   *
   */ 
  std::function<double(double)> beta_x = 
    [] (double x)-> double{
      (void)x; 
      return -r*x; 
    };

  /**
   * @function initial_cond
   *
   * @brief Evaluate the initial condition at 'x', 
            scalar-valued
   *
   * @param x the coordinates of the point to evaluate the function at
   * @return a scalar the value of the initial condition at 'x'
   *
   */ 
  std::function<double(double)> initial_cond = 
   [](double x)->double{
     if (x>10)
       return 1;
     else
       return 0;
   };

  /**
   * @function f_x
   *
   * @brief Evaluate the load functor f(x), 
            scalar-valued
   *
   * @param x the coordinates of the point to evaluate the function at
   * @return a scalar the value of f at 'x'
   *
   */ 
  std::function<double(double x)> f_x =
    [] (double x) -> double{
      (void)x; 
      return 0;
    };

} // namespace closure

#endif
