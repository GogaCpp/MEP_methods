
#include <iostream>
#include <fstream>
#include "mepClass.cpp"
#include <numbers>

double DWcircle(const std::vector<double> arr)
{
	double x=arr[0];
	double y = arr[1];
	return (1-x*x-y*y)*(1-x*x-y*y) + y*y /(x*x+y*y);
}

double SC(const std::vector<double>& x)
{
  double a=x[0],b=x[1];
  double l =std::sin(1.3*a)*std::cos(.9*b)+std::cos(.8*a)*std::sin(1.9*b)+std::cos(a*.2*b);
  return l;
  }

double MullerBrowns(const std::vector<double> x)
{
	double A[4] = {-200, -100, -170, 15};
	double a[4] = {-1, -1, -6.5, 0.7};
	double b[4] = {0, 0, 11, 0.6};
	double c[4] = {-10, -10, -6.5, 0.7};
	double x0[4] = {1, 0, -0.5, -1};
	double y0[4] = {0, 0.5, 1.5, 1};
	double res = 0;
	for(int i = 0; i < 4; i++)
	{
		res += A[i] * std::exp(a[i] * (x[0] - x0[i]) * (x[0] - x0[i]) + b[i] * ((x[0] - x0[i]) * (x[1] - y0[i])) + c[i] * (x[1] - y0[i]) * (x[1] - y0[i]));
	}
	return res;
}

double DWsimple(const std::vector<double> arr)
{
	double x=arr[0];
	double y = arr[1];
	return (x*x-1)*(x*x-1) + y*y;
}
double  Schwefel(const std::vector<double> arr)
{
  double result=0;
  for(int i=0;i<arr.size();i++)
    result = std::max(std::abs(arr[i]), result);
  return result;
  }
int main()
{
  

	
  std::vector<double> head = {0, 0, 0}; 
  std::vector<double> tail = {-4.1 , -4.2, -4.3};
  mep::STRING<double> StringMethod;
	StringMethod.set_parameters(head,tail,60,0.01,10000,1.0e-6,1.0);
	StringMethod.set_flexibility(1);
	StringMethod.set_function(DWcircle);
	//StringMethod.set_optimal_path(op);
  StringMethod.compute_mep();
	std::vector<std::vector<double>> path_str = StringMethod.get_mep();
  print2file2d("string.dat",path_str);
  
  mep::NEB<double> NebMethod;
	NebMethod.set_parameters(head,tail,60,0.01,10000,1.0e-6,1.0);
	NebMethod.set_flexibility(1);
	NebMethod.set_function(DWcircle);
	//NebMethod.set_optimal_path(op);
  NebMethod.compute_mep();
	std::vector<std::vector<double>> path_neb = NebMethod.get_mep();
  print2file2d("neb.dat",path_neb);
  
 
  }
