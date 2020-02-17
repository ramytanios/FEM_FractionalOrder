#ifndef ROUTINES_HPP
#define ROUTINES_HPP 

#include <fstream>
#include <math.h>

/*
 * @struct log_op
 *
 * @brief A wrapper around the std::log function,
 *	  mainly used for std::transform, std::for_each
 *
 * @param arg a double to evaluate the logarithm at
 * @return the logarithm of of arg
 *
 */
struct log_op{
  double operator()(double arg){
    return std::log(arg);
  }
};

// template<class _MatrixType>
// inline double conditionNumber(_MatrixType& A){
//   Eigen::JacobiSVD<_MatrixType> svd(A); 
//   singularValues= svd.singularValues(); 
//   return singularValues(0) / 
//     singularValues(singularValues.size()-1); 
// }
// 
// template<class _MatrixType> 
// std::vector<double> conditionNumberSinc(_MatrixType& L, _MatrixType& A,
//     		double h, double beta){





/**
 * @function getIndicesOfInterest
 *
 * @brief Returning the indices of a vector  for a given domain
 *	  of interest, taking advantage of the symmetry of intervals
 *
 * Computational domain = [-R,R] ([-4,4]) 
 * Domain of interest = [lowerDomOfInt, upperDomOfInt] ([-2,2])
 *
 * @throws std::invalid_argument thrown if 'lowerDomOfInt' is positive
 * @param x the vector of interest
 * @param lowerDomOfInt the lower bound of the domain of interest
 * @return the indices corresponding to the bounds of the domain of interest
 *
 */
std::pair<size_t, size_t> getIndicesOfInterest(Eigen::VectorXd& x,
    		int lowerDomOfInt){ 
  if (lowerDomOfInt >= 0) 
    throw std::invalid_argument("lowerDomOfInt must be negative"); 
  int i;
  for (i=0; i<x.size(); i++){ 
    if(x[i] >= lowerDomOfInt)
      break;
  }
  return {i, x.size()-1-i}; 
}
/* */

/**
 * @function ReadNthLine
 *
 * @brief Reading the n-th line from a text file
 *
 * @throws std::runtime_error if cannot read the text file
 * @param filename the name of the file to read from
 * @param N the number of the line to read
 * @param size the size of the vector where to store the read line
 * @return the vector containing the read line
 *
 */
Eigen::VectorXd ReadNthLine(const std::string& filename, int N, size_t size){
  std::ifstream in; 
  in.exceptions(std::ifstream::failbit); 
  try{ 
    in.open(filename);
  }
  catch (const std::exception& e){ 
    std::ostringstream msg;
    msg << "File " << filename << " cannot be opened. "; 
    throw std::runtime_error(msg.str());
  }
  std::string s;
  s.reserve(300); /// for performance
  for(int i = 0; i < N; ++i) /// skip first N-1 lines 
      std::getline(in, s);
  std::getline(in,s); /// read as a string
  Eigen::VectorXd res(size); ;
  std::stringstream ss(s); 
  for (int i=0; i<size; i++)
    ss >> res[i]; 
  return res; 
}
/*  www.stackoverflow.com */ 

/**
 * @function linear_fit
 *
 * @brief Performing 1D least squares linear regression, to approximate
 * 	  the algebraic rate of convergence
 *
 * @tparam vectorType the type of the input vectors
 * @param y_std,x_std the vectors to be fitted against each other
 * @return the coefficients of the linear regression, intercept and slope.
 *
 */
template<typename vectorType> 
std::pair<double,double> linear_fit(vectorType& y_std, vectorType& x_std){
  if (x_std.size() != y_std.size())
    throw std::length_error("The two input vectors should be of same length!"); 
 
  Eigen::VectorXd y = 
    Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(y_std.data(), y_std.size()); 
  Eigen::VectorXd x = 
    Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(x_std.data(), x_std.size()); 
  Eigen::MatrixXd X(x.size(), 2); 
 
  X << Eigen::VectorXd::Constant(x.size(), 1), x; 
  Eigen::VectorXd sol = X.fullPivHouseholderQr().solve(y);
 
  return std::make_pair(sol[0], sol[1]); 
}
/* */ 

/** 
 * @function l2_norm 
 *
 * @brief Calculating the L2-norm of a one dimensional 
 * 	  FEM function using the local trapezoidal rule. 
 *
 * @tparam arg the coefficients basis expansion of the function
 * @tparam mesh_size the mesh size
 * @return the L2 norm of the function described by arg
 *
 */ 
double l2_norm(Eigen::VectorXd& arg, double mesh_size){ 
  double result; 
  size_t length = arg.size(); 
  result = std::pow(arg(0),2) + std::pow(arg(length-1),2); 
  for (int i=1; i<length-1; i++){ 
    result += 2*std::pow(arg(i),2); ;
  }
  result = sqrt(mesh_size/2.0 * result); 
  return result;
}
/* */ 

/** 
 * @function N
 *
 * @brief Evaluates the CDF of a standard normal
 *	  at given points. 
 *
 * @tparam arg vector of values where to evaluate the CDF
 * @return a vector of evaluations
 *
 */ 
Eigen::VectorXd N(Eigen::VectorXd& arg){
  Eigen::VectorXd result(arg.size()); 
  std::transform(arg.data(), arg.data()+arg.size(), 
      result.data(), [](double s){
      return 0.5 * (1+erf(s/std::sqrt(2)));}); 
  return result; 
}

/** 
 * @function CallBS 
 *
 * @brief Calculates the price of a call option
 *	  according to Black-Scholes formula.
 *
 * @tparam S a vector of prices
 * @tparam K the strike price
 * @tparam r the risk-free interest rate
 * @tparam sigma the volatility
 * @tparam t the time t
 * @tparam T the maturity time
 * @return a vector of prices according to BS formula
 *
 */ 

Eigen::VectorXd CallBS(Eigen::VectorXd& S, double K, double r,
    double sigma, double t, double T){
  
  Eigen::VectorXd price(S.size()); 
  double tau = T-t; 
  
  Eigen::VectorXd d1 = ( S.array().log() - std::log(K) + (r+0.5 * std::pow(sigma,2)) * tau ) 
    /sigma/ std::sqrt(tau); 
  Eigen::VectorXd d2 = d1.array() - sigma*std::sqrt(tau); 

  price = S.array() * N(d1).array() - K*std::exp(-r*tau)*N(d2).array(); 

  return price; 
}

/** 
 * @function PutBS 
 *
 * @brief Calculates the price of a put option
 *	  according to Black-Scholes formula.
 *
 * @tparam S a vector of prices
 * @tparam K the strike price
 * @tparam r the risk-free interest rate
 * @tparam sigma the volatility
 * @tparam t the time t
 * @tparam T the maturity time
 * @return a vector of prices according to BS formula
 *
 */ 


Eigen::VectorXd PutBS(Eigen::VectorXd& S, double K, double r,
    double sigma, double t, double T){
  
  Eigen::VectorXd price(S.size()); 
  double tau = T-t; 
  
  Eigen::VectorXd d1 = -( S.array().log() - std::log(K) + (r+0.5 * std::pow(sigma,2)) * tau ) 
    /sigma/ std::sqrt(tau); 
  Eigen::VectorXd d2 = d1.array() + sigma*std::sqrt(tau); 

  price = -S.array() * N(d1).array() + K*std::exp(-r*tau)*N(d2).array(); 

  return price; 
}

/** 
 * @function DigitalBS 
 *
 * @brief Calculates the price of a digital option
 *	  according to Black-Scholes formula.
 *
 * @tparam S a vector of prices
 * @tparam K the strike price
 * @tparam r the risk-free interest rate
 * @tparam sigma the volatility
 * @tparam t the time t
 * @tparam T the maturity time
 * @return a vector of prices according to BS formula
 *
 */ 
Eigen::VectorXd DigitalBS(Eigen::VectorXd& S, double K, double r,
    double sigma, double t, double T){
  
  Eigen::VectorXd price(S.size()); 
  double tau = T-t; 
  
  Eigen::VectorXd d1 = ( -S.array().log() + std::log(K) - (r-0.5 * std::pow(sigma,2)) * tau ) 
    /sigma/ std::sqrt(tau); 

  price = std::exp(-r*tau)*(1-N(d1).array()); 

  return price; 
}

#endif 
