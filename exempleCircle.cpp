
#include <iostream>
#include <fstream>
#include "mepClass.cpp"


double DWcircle(const std::vector<double> arr)
{
	double x=arr[0];
	double y = arr[1];
	return (1-x*x-y*y)*(1-x*x-y*y) + y*y /(x*x+y*y);
}

int main()
{
  
  std::vector<std::vector<double>> op={
    {-1.000000e+00, 8.979318e-11},
    {-9.396926e-01, 3.420201e-01},
    {-7.660444e-01, 6.427876e-01},
    {-5.000000e-01, 8.660254e-01},
    {-1.736482e-01, 9.848078e-01},
    {1.736482e-01,  9.848078e-01},
    {5.000000e-01,  8.660254e-01},
    {7.660444e-01,  6.427876e-01},
    {9.396926e-01,  3.420201e-01},
    {1.000000e+00,  0.000000e+00}
    };
	std::vector<double> head = {-0.5, 0.5}; 
  std::vector<double> tail = {0.5 , 0.5};
  mep::STRING StringMethod;
	StringMethod.set_parameters(head,tail);
	StringMethod.set_flexibility(0);
	StringMethod.set_function(DWcircle);
	StringMethod.set_optimal_path(op);
  StringMethod.compute_mep();
	std::vector<std::vector<double>> path_str = StringMethod.get_mep();
  print2file2d("string_circle_10.dat",path_str);
  
  mep::NEB NebMethod;
	NebMethod.set_parameters(head,tail);
	NebMethod.set_flexibility(0);
	NebMethod.set_function(DWcircle);
	NebMethod.set_optimal_path(op);
  NebMethod.compute_mep();
	std::vector<std::vector<double>> path_neb = NebMethod.get_mep();
  print2file2d("neb_circle_10.dat",path_neb);
  
  }
