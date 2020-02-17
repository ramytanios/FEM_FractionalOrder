#ifndef PDECOEFF_HPP
#define PDECOEFF_HPP

namespace PdeCoefficients{

  /**
   * @function alpha_x
   *
   * @brief Evaluate the coefficient functor alpha(x), 
            matrix-valued
   *
   * @param x the coordinates of the point to evaluate the function at
   * @return a 2x2 matrix the value of alpha at 'x'
   *
   */ 
  std::function<Eigen::Matrix2d(Eigen::Vector2d)> alpha_x =
    [] (Eigen::Vector2d x)->Eigen::Matrix2d{
        Eigen::Matrix2d res;
	(void)x; 
	  res << 1,0,
        	 0,1;
	  return res; 
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
  std::function<double(Eigen::Vector2d)> gamma_x = 
    [] (Eigen::Vector2d x)-> double{
      (void)x; 
      return 0;
    };
  
  /**
   * @function beta_x
   *
   * @brief Evaluate the coefficient functor beta(x), 
            vector-valued
   *
   * @param x the coordinates of the point to evaluate the function at
   * @return a 2x1 vector the value of beta at 'x'
   *
   */ 
  std::function<Eigen::Vector2d(Eigen::Vector2d)> beta_x = 
    [] (Eigen::Vector2d x)-> Eigen::Vector2d{
      (void)x; 
      Eigen::Vector2d result; 
      result << 0,0; 
      return result;
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
  std::function<double(Eigen::Vector2d)> initial_cond = 
   [](Eigen::Vector2d x)->double{
     return x[0]*(1-x[0])*x[1]*(1-x[1]); 
//     return 5; 
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
  std::function<double(Eigen::Vector2d x)> f_x =
    [] (Eigen::Vector2d x) -> double{
      (void)x; 
      return 0;
    };

} // namespace closure

#endif
