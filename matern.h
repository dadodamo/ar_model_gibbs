#ifndef MATERN_H
#define MATERN_H

#include <iostream>
#include<cmath>
#include <boost/math/special_functions/bessel.hpp>


// distance function and Gamma function definition
// bessel function from boost

static double matern(double& x, double& phi, double& nu){
    double ans;
    if(x == 0)
       ans = 1;
    else
       ans = 1.0 / ( exp(std::lgamma(nu)) * pow(2, nu - 1.0) ) * pow(x, nu) * boost::math::cyl_bessel_k(nu, x);
    return ans;
 }
// LOG scale
#endif