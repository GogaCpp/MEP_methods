
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
	//std::vector<double> head = {-0.5, 0.5}; 
  //std::vector<double> tail = {0.5 , 0.5};
  std::vector<double> head = {0, 0, 0}; 
  std::vector<double> tail = {-4.1 , -4.2, -4.3};
  mep::STRING<double> StringMethod;
	StringMethod.set_parameters(head,tail,60,0.01,10000,1.0e-6,1.0);
	StringMethod.set_flexibility(1);
	StringMethod.set_function(SC);
	//StringMethod.set_optimal_path(op);
  StringMethod.compute_mep();
	std::vector<std::vector<double>> path_str = StringMethod.get_mep();
  print2file2d("string_rs_coef.dat",path_str);
  //print2file1d("string_sw_rmse.dat",StringMethod.get_rmse_evolution());
  mep::NEB<double> NebMethod;
	NebMethod.set_parameters(head,tail,60,0.01,10000,1.0e-6,1.0);
	NebMethod.set_flexibility(1);
	NebMethod.set_function(SC);
	//NebMethod.set_optimal_path(op);
  NebMethod.compute_mep();
	std::vector<std::vector<double>> path_neb = NebMethod.get_mep();
  print2file2d("neb_rs.dat",path_neb);
  //print2file1d("neb_sw_rmse.dat",NebMethod.get_rmse_evolution());
  
  /*mep::NEB<double> NebMethod1;
	NebMethod1.set_parameters(head,tail,25,0.0001,10000,1.0e-6,0.5);
	NebMethod1.set_flexibility(1);
	NebMethod1.set_function(rastrigin_simple);
	//NebMethod.set_optimal_path(op);
  NebMethod1.compute_mep();
	std::vector<std::vector<double>> path_neb1 = NebMethod1.get_mep();
  print2file2d("INIT_SC.dat",NebMethod1.get_initial_path());
  
  print2file1d("neb_SC_rmse05.dat",NebMethod1.get_rmse_evolution());*/
  
  
  
  /*std::vector<double> head = {0, 0}; 
  std::vector<double> tail = {-1, 0};
  for(int i=0;i<50;i++)
	{
	  mep::NEB<double> obj;
	  obj.set_parameters(head,tail,10,0.01,10000);
	  obj.set_flexibility(0);
	  obj.set_function(DWcircle);
	  //obj.set_optimal_path(op);
	  auto time_one = std::chrono::steady_clock::now();
	  std::vector<std::vector<double>> pa = obj.compute_mep();
    auto time_two = std::chrono::steady_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(time_two - time_one).count();
    print2file2d("timeCS.dat",pa);
    std::cout << elapsed_time<<", " ;

	  }*/
  }
