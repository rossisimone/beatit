/*
 * CTestUtil.cpp
 *
 *  Created on: Sep 6, 2016
 *      Author: srossi
 */

#include "Util/CTestUtil.hpp"
#include <iostream>
#include <cmath>

namespace BeatIt
{

namespace CTest
{


int check_test(double norm, double reference_norm, double tol)
{
    double diff = std::abs(norm-reference_norm);
    if(diff < tol)
	{
  	  std::cout << "Well done! Test was succesful! Diff = " << diff <<  std::endl;
  	  return EXIT_SUCCESS;
	}
    else
	  {
    	  std::cout << "Failure: diff  = " << diff  << std::endl;
    	  std::cout << "Norm is  = " << norm << ", but was "  << reference_norm  << std::endl;
    	  return EXIT_FAILURE;
	  }
}

} //CTest

} //BeatIt
