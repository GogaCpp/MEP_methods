#include "amep.h"
#include <iostream>
#include <chrono>
#include <fstream>

void printString(const std::vector<std::vector<double>>& string);

namespace mep
{
	class NEB : IMEP
	{
	private:
		vect head;
		vect tail;
		vect rmse_evolution;
		vect profile;
		path MEP;
		path optimal_path;
		path initial_path;
		std::vector<path> mep_evolution;
		size_t point_number;
		size_t flex;
		size_t niter;
		double step;
		double h;
		double k;
		std::function<double(const vect&)> mainfunc;


	public:
		void set_parameters(const vect& _head,
							const vect& _tail,
							size_t pn = 8,
							double step = 0.01,
							size_t niter = 2500,
							double h = 1.0e-6,
							double k = 1.0)
		{
			this->head=_head;
			this->tail=_tail;
			this->point_number=pn;
			this->step=step;
			this->niter=niter;
			this->h=h;
			this->k=k;
		}
		void set_step(double dt=1.0e-3)
		{
			this->step=dt;
		}

		void set_function(std::function<double(const vect&)> fn)
		{
			mainfunc = fn;
		}
		void set_flexibility(size_t flex)
		{
			this->flex=flex;
		}
		path get_mep()
		{
			return MEP;
		}
		void set_optimal_path(const path& optimal_path)
		{
			this->optimal_path=optimal_path;
		}
		std::vector<double> get_profile() const
		{
			return profile;

		}
		double get_rmse(const path& a, const path& b) const
		{
			if(a.size() != b.size())
			{
				std::cout<<"Вектора для подсчета RMSE разных размеров"<<std::endl;
				return -1;
			}
			else
			{
				double sum = 0.0;
				for(size_t i = 0; i < a.size(); i++)
				{
					for(size_t j = 0; j < a[i].size(); j++)
					{
						double diff = a[i][j] - b[i][j];
						sum += diff * diff;
					}
				}
				double mse = sum / (a.size()*a.front().size());
				return std::sqrt(mse);
			}
		}
		void set_initial_path(const path& a)
		{
			this->initial_path = a;
		}

		vect get_rmse_evolution()
		{
			return rmse_evolution;
		}
		std::vector<path> get_mep_evolution()
		{
			return mep_evolution;
		}
		path get_initial_path() const
		{
			return this->initial_path;
		}

		void numericalGradient(double x, double y, double& grad_x, double& grad_y)
		{
			double h = 1.0e-6;
			double fx = mainfunc({x, y});
			double fxh1 = mainfunc({x + h, y});
			double fxh2 = mainfunc({x, y + h});

			grad_x = (fxh1 - fx) / h;  // Приближенная производная по x
			grad_y = (fxh2 - fx) / h;  // Приближенная производная по y
		}
		double score(const std::vector<double>& v)
		{
			//return -exp(-((v[0]-1)*(v[0]-1) + (v[1]-1)*(v[1]-1))) - exp(-((v[0]+1)*(v[0]+1) + (v[1]+1)*(v[1]+1))); //MullerBrowns(v[0], v[1]);
			return mainfunc({v[0], v[1]});
		}
		double grad(const std::vector<double> &vct, size_t n)
		{
			double h = 1.0e-6;
			std::vector<double> vp = vct;
			std::vector<double> vm = vct;
			vp[n] += h;
			vm[n] -= h;
			return	(score(vp) - score(vm)) / (2*h);
		}
		vect findLocalMinimum(vect x, double alpha=0.000001, int maxIterations= 100000)
		{
			vect grad_xyz(x.size());
			for(int i = 0; i < maxIterations; ++i)
			{

				for(int j=0; j<x.size(); j++)
				{
					grad_xyz[j]=get_gradient(x,j);
					x[j]-=alpha*grad_xyz[j];
				}
			}
			return x;
		}
		double scalar_product(const std::vector<double>& vector1, const std::vector<double>& vector2)
		{
			double result = 0.0;
			for(size_t i = 0; i < vector1.size(); ++i)
				result += vector1[i] * vector2[i];
			return result;
		}
		void initialize_line(std::vector<std::vector<double>>& str, const std::vector<double>& head, const std::vector<double>& tail)
		{
			std::vector<double> es(head.size());
			for(size_t i = 0; i < es.size(); i++)
				es[i] = tail[i] - head[i];
			size_t n = str.size() - 1;
			for(size_t i = 0; i < str.size(); i++)
			{
				for(size_t j = 0; j < str[i].size(); j++)
				{
					str[i][j] = head[j] + i*es[j]/n;
				}
			}
		}
		void get_tangent(std::vector<std::vector<double>>& images, const std::vector<std::vector<double>>& str)
		{
			for(size_t i = 1; i != images.size() - 1; i++)
			{

				std::vector<double> values(images[i].size());
				double sum = 0.0;
				for(size_t j = 0; j < images[i].size(); j++)
				{
					values[j] = str[i+1][j] - str[i-1][j];
					sum += values[j]*values[j];
				}
				double norm = std::sqrt(sum);
				for(size_t j = 0; j < images[i].size(); j++)
				{
					images[i][j] = values[j]/norm;
				}
			}
		}

		void get_parallel(std::vector<std::vector<double>>& images, const std::vector<std::vector<double>>& str,const std::vector<std::vector<double>>& tangent, double k)
		{
			std::vector<std::vector<double>> amp(images.size(), std::vector<double>(2));
			for(size_t i = 1; i != images.size() - 1; i++)
			{
				double sum = 0.0;
				for(size_t j = 0; j < images[i].size(); j++)
				{
					sum += (str[i-1][j] - str[i][j])*(str[i-1][j] - str[i][j]);
				}
				amp[i][0] = std::sqrt(sum);

				sum = 0.0;
				for(size_t j = 0; j < images[i].size(); j++)
				{
					sum += (str[i+1][j] - str[i][j])*(str[i+1][j] - str[i][j]);
				}
				amp[i][1] = std::sqrt(sum);
			}
			std::cout.precision(16);
			for(size_t i = 0; i < images.size(); i++)
			{
				double diff = amp[i][1] - amp[i][0];
				for(size_t j = 0; j < images[i].size(); j++)
				{
					images[i][j] = k * diff * tangent[i][j];
				}
			}
		}
		double get_gradient(const vect& v, size_t n,double h = 0.0000001) const
		{
			vect a = v;
			vect b = v;
			a[n]+=h;
			b[n]-=h;
			return (mainfunc(a) - mainfunc(b))/(2*h);
		}

		void get_gradients(std::vector<std::vector<double>>& images, const std::vector<std::vector<double>>& str)
		{
			for(size_t i = 1; i < images.size()-1; i++)
			{
				for(size_t j = 0; j < images[j].size(); j++)
				{
					images[i][j] = grad(str[i], j);
				}
			}
		}
		void get_prepared(std::vector<std::vector<double>>& images, const std::vector<std::vector<double>>& tangent)
		{
			for(size_t i = 1; i < images.size()-1; i++)
			{
				double dot = scalar_product(images[i], tangent[i]);
				for(size_t j = 0; j < images[j].size(); j++)
				{
					images[i][j] = images[i][j] - dot*tangent[i][j];
				}
			}
		}
		void compute_mep()
		{
			path images(point_number+2, std::vector<double>(head.size()));

			for(size_t i=0; i<head.size(); i++)
			{
				images[0][i]=head[i];
				images[images.size()-1][i] = tail[i];
			}

			std::vector<std::vector<double>> tangents(point_number+2, std::vector<double>(head.size()));
			std::vector<std::vector<double>> springs(point_number+2, std::vector<double>(head.size()));
			std::vector<std::vector<double>> gradients(point_number+2, std::vector<double>(head.size()));

			std::vector<double> force(point_number*head.size(), 0);
			std::vector<double> old_force(point_number*head.size(), 0);
			std::vector<double> velocity(point_number*head.size(), 0);
			if(initial_path.empty())
			{
				initialize_line(images, head, tail);
				initial_path = images;
			}
			else
				images = initial_path;
			mep_evolution.push_back(images);

			for(size_t g = 0; g < niter; g++)
			{
				get_tangent(tangents, images);
				get_parallel(springs, images, tangents, k);
				get_gradients(gradients, images);
				get_prepared(gradients, tangents);


				for(size_t i = 1, l = 0; i < springs.size() - 1; i++)
				{
					for(size_t j = 0; j < springs[i].size(); j++, l++)
					{
						force[l] = springs[i][j] - gradients[i][j];
					}
				}

				auto vf = scalar_product(velocity, force);
				if(vf > 0.0)
				{
					double sum = 0.0;
					for(size_t i = 0; i < force.size(); i++)
						sum += force[i]*force[i];
					for(size_t i = 0; i < force.size(); i++)
						velocity[i] = vf*force[i]/sum;
				}
				else
				{
					for(size_t i = 0; i < velocity.size(); i++)
						velocity[i] = 0.0;
				}

				for(size_t i = 0; i < velocity.size(); i++)
				{
					velocity[i] += step*(old_force[i] + force[i])/2.0;
				}

				for(size_t i = 1; i < images.size()-1; i++)
				{
					for(size_t j = 0; j < images[i].size(); j++)
					{
						images[i][j] += velocity[head.size()*(i-1) + j]*step;
					}
				}
				if(this->flex == 0)
				{
					images[0] = findLocalMinimum(images[0],h,200);
					images[images.size()-1]=findLocalMinimum(images[images.size()-1],h,200);
				}
				else if(this->flex == 2)
				{
					images[images.size()-1]=findLocalMinimum(images[images.size()-1],h,200);

				}
				else if(this->flex == 3)
				{
					images[0] = findLocalMinimum(images[0],h,200);
				}


				for(size_t i = 0; i < old_force.size(); i++)
				{
					old_force[i] = force[i];
				}
				mep_evolution.push_back(images);
				if(!optimal_path.empty())
				{
					rmse_evolution.push_back(this->get_rmse(optimal_path,images));
				}
				else
				{
					rmse_evolution.push_back(this->get_rmse(mep_evolution[g],images));
				}
			}
			this->MEP= images;
		}
	};
	class STRING : IMEP
	{
	private:
		vect head;
		vect tail;
		vect rmse_evolution;
		vect profile;
		path MEP;
		path optimal_path;
		path initial_path;
		std::vector<path> mep_evolution;
		size_t point_number;
		size_t flex;
		size_t niter;
		double step;
		double h;
		double k;
		std::function<double(const vect&)> mainfunc;


	public:
		void set_parameters(const vect& _head,
							const vect& _tail,
							size_t pn = 10,
							double step = 0.01,
							size_t niter = 10000,
							double h = 1.0e-6,
							double k = 1.0)
		{

			this->head=_head;
			this->tail=_tail;
			this->point_number=pn;
			this->step=step;
			this->niter=niter;
			this->h=h;
			this->k=k;
		}
		void set_function(std::function<double(const vect&)> fn)
		{
			mainfunc = fn;
		}
		void set_flexibility(size_t flex)
		{
			this->flex=flex;
		}
		path get_mep()
		{
			return MEP;
		}
		void set_optimal_path(const path& optimal_path)
		{
			this->optimal_path=optimal_path;
		}
		std::vector<double> get_profile() const
		{
			return profile;
		}
		double get_rmse(const path& a, const path& b) const
		{
			if(a.size() != b.size())
			{
				/*std::cout<<"Вектора для подсчета RMSE разных размеров"<<std::endl;*/
				return -1;
			}
			else
			{
				double sum = 0.0;
				for(size_t i = 0; i < a.size(); i++)
				{
					for(size_t j = 0; j < a[i].size(); j++)
					{
						double diff = a[i][j] - b[i][j];
						sum += diff * diff;
					}
				}
				double mse = sum / (a.size()*a.front().size());
				return std::sqrt(mse);
			}
		}
		void set_initial_path(const path& a)
		{
			this->initial_path = a;
		}
		vect get_rmse_evolution()
		{
			return rmse_evolution;
		}
		std::vector<path> get_mep_evolution()
		{
			return mep_evolution;
		}
		path get_initial_path() const
		{
			return this->initial_path;
		}

		void initialize_line(std::vector<std::vector<double>>& str, const std::vector<double>& head, const std::vector<double>& tail)
		{
			std::vector<double> es(head.size());
			for(size_t i = 0; i < es.size(); i++)
				es[i] = tail[i] - head[i];
			size_t n = str.size() - 1;
			for(size_t i = 0; i < str.size(); i++)
			{
				for(size_t j = 0; j < str[i].size(); j++)
				{
					str[i][j] = head[j] + i*es[j]/n;
				}
			}
		}
		std::vector<double> linear_interpolation(const std::vector<double>& gamma,
				const std::vector<double>& xp,
				const std::vector<double>& fp)
		{
			std::vector<double> interpolated_values;

			for(size_t i = 0; i < gamma.size(); i++)
			{
				auto lower = std::lower_bound(xp.begin(), xp.end(), gamma[i]);
				size_t idx = std::max(static_cast<int>(lower - xp.begin()) - 1, 0);

				double x1 = xp[idx];
				double x2 = xp[idx + 1];
				double y1 = fp[idx];
				double y2 = fp[idx + 1];

				double interpolated_y = y1 + (gamma[i] - x1) * (y2 - y1) / (x2 - x1);
				interpolated_values.push_back(interpolated_y);
			}

			return interpolated_values;
		}

		double norm(const std::vector<double>& v)
		{
			double sum = 0.0;
			for(double val : v)
			{
				sum += val * val;
			}
			return std::sqrt(sum);
		}


		double get_gradient(const vect& v, size_t n,double h = 0.0000001) const
		{
			vect a = v;
			vect b = v;
			a[n]+=h;
			b[n]-=h;
			return (mainfunc(a) - mainfunc(b))/(2*h);
		}

		std::vector<std::vector<double>> step_euler(std::vector<std::vector<double>>& string, double step)
		{

			path string_grad(string.size(),vect(string[0].size()));
			for(size_t i = 0; i < string.size(); i++)
			{
				for(size_t j=0; j< string[0].size(); j++)
					string_grad[i][j]=get_gradient(string[i],j);

			}
			double h = 0.0;
			for(size_t i = 0; i < string_grad.size(); i++)
			{
				double res=0;
				for(size_t j=0; j< string_grad[j].size(); j++)
					res+=string_grad[i][j]*string_grad[i][j];
				h = std::max(h, std::sqrt(res));

			}

			size_t i = 0;
			size_t en = string.size();

			if(flex==1) {i++; en--;}
			if(flex==2) {i++;}
			if(flex==3) {en--;}

			for(i; i<en; i++)
			{
				for(size_t j=0; j<string[0].size(); j++)
					string[i][j]-=step * string_grad[i][j] / h;
			}
			return string;
		}

		std::vector<double> linspace(double start, double end, int npts)
		{
			std::vector<double> result;
			double step = (end - start) / (npts - 1);
			for(int i = 0; i < npts; ++i)
			{
				result.push_back(start + i * step);
			}
			return result;
		}

		void compute_mep()
		{

			path string_xyz(point_number,std::vector<double>(head.size()));
			initialize_line(string_xyz,head,tail);
			std::vector<std::vector<double>> string;
			if(initial_path.empty()) { string=string_xyz; initial_path=string_xyz;}
			else   string = initial_path;
			mep_evolution.push_back(string);

			std::vector<std::vector<double>> old_string(string.size(), std::vector<double>(string[0].size(), 0.0));


			for(int i=1; i<niter+1; i++)
			{

				old_string=string;

				string=step_euler(string,step);

				std::vector<double> arclength;
				arclength.push_back(0.0);

				for(size_t it = 1; it < string.size(); it++)
				{
					std::vector<double> diff;
					for(size_t j = 0; j < string[it].size(); j++)
					{
						diff.push_back(string[it][j] - string[it - 1][j]);
					}
					double norm_diff = norm(diff);
					arclength.push_back(arclength.back() + norm_diff);
				}
				double arclength_last = arclength.back();
				for(double &value : arclength)
				{
					value /= arclength_last;
				}

				path xp_xyz(string[0].size(),std::vector<double>(arclength));
				std::vector<double> gamma=linspace(0,1,point_number);
				path fp_xyz(string[0].size(),std::vector<double>(arclength));
				path reparam_xyz(string[0].size(),std::vector<double>(string.size()));


				for(size_t itera =0 ; itera<string[0].size(); itera++)
				{
					xp_xyz[itera]=arclength;
					for(size_t k=0; k<string.size(); k++)
					{
						fp_xyz[itera][k]=string[k][itera];
					}
				}

				for(size_t itera =0 ; itera<string[0].size(); itera++)
				{

					reparam_xyz[itera] =linear_interpolation(gamma, xp_xyz[itera], fp_xyz[itera]);
				}
				if(string.size() != reparam_xyz[0].size() || string[0].size() != reparam_xyz.size())
				{
					std::cout << "Error: Sizes of string and reparam vectors do not match." << std::endl;
					return ;
				}


				for(size_t k = 0; k < string.size(); k++)
				{
					for(size_t it=0; it<string[0].size(); it++)
						string[k][it]=reparam_xyz[it][k];

				}
				this->mep_evolution.push_back(string);

				if(!optimal_path.empty())
				{
					this->rmse_evolution.push_back(this->get_rmse(optimal_path,string));
				}
				else
				{
					this->rmse_evolution.push_back(this->get_rmse(mep_evolution[i-1],string));
				}
			}
			this->MEP = string;
			profile = vect(MEP.size());
			for(size_t i =0; i<this->MEP.size(); i++)
				this->profile[i]=(mainfunc(this->MEP[i]));
		}

	};
}

void print2file1d(std::string fname, const std::vector<double> &u, size_t step = 1)
{
	std::ofstream file(fname);
	if(file.is_open())
	{
		for(size_t i = 0; i < u.size(); i+=step)
			file <<std::scientific<<i <<'\t' <<u[i]<<std::endl;
		//std::cout<<"Вывод в фаил закончен"<<std::endl;
		file.close();
	}
}

void print2file2d(std::string fname, const std::vector<std::vector<double>> &u, size_t step = 1)
{
	std::ofstream file(fname);
	if(file.is_open())
	{
		for(size_t i = 0; i < u.size(); i+=step)
		{
			for(size_t j=0; j<u[i].size(); j++)
				file <<std::scientific<< u[i][j]<<'\t';
			file<<std::endl;
		}
		//std::cout<<"Вывод в фаил закончен"<<std::endl;
		file.close();
	}
}


