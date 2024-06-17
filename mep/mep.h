#ifndef MEP_H
#define MEP_H

#include <amep.h>
#include <iostream>
#include <chrono>
#include <fstream>
void printString(const std::vector<std::vector<double>>& string);

namespace mep
{
  template <typename U>
	class NEB :  public IMEP<U>
	{
    using path = std::vector<std::vector<U>>;
		using vect = std::vector<U>;
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
		U step;
		U h;
		U k;
		std::function<U(const vect&)> mainfunc;


	public:
		void set_parameters(const vect& _head,
							const vect& _tail,
							size_t pn = 10,
							U step = 0.01,
							size_t niter = 2500,
							U h = 1.0e-6,
							U k = 1.0)
		{
			this->head=_head;
			this->tail=_tail;
			this->point_number=pn-2;
			this->step=step;
			this->niter=niter;
			this->h=h;
			this->k=k;
		}
		

		void set_function(std::function<U(const vect&)> fn)
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
		std::vector<U> get_profile() const
		{
			return profile;

		}
		U get_rmse(const path& a, const path& b) const
		{
				U sum = 0.0;
				for(size_t i = 0; i < a.size(); i++)
				{
					for(size_t j = 0; j < a[i].size(); j++)
					{
						U diff = a[i][j] - b[i][j];
						sum += diff * diff;
					}
				}
				U mse = sum / (a.size()*a.front().size());
				return std::sqrt(mse);
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

	
		vect findLocalMinimum(vect x, U alpha=0.000001, int maxIterations= 100000)
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
		U scalar_product(const std::vector<U>& vector1, const std::vector<U>& vector2)
		{
			U result = 0.0;
			for(size_t i = 0; i < vector1.size(); ++i)
				result += vector1[i] * vector2[i];
			return result;
		}
    U norm(const std::vector<U>& v)
    {
      return std::sqrt(scalar_product(v,v));
      }
		void initialize_line(std::vector<std::vector<U>>& str, const std::vector<U>& head, const std::vector<U>& tail)
		{
			std::vector<U> es(head.size());
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
		void get_tangent(std::vector<std::vector<U>>& images, const std::vector<std::vector<U>>& str)
		{
			for(size_t i = 1; i != images.size() - 1; i++)
			{

				std::vector<U> values(images[i].size());
				U sum = 0.0;
				for(size_t j = 0; j < images[i].size(); j++)
				{
					values[j] = str[i+1][j] - str[i-1][j];
					sum += values[j]*values[j];
				}
				U norm = std::sqrt(sum);
				for(size_t j = 0; j < images[i].size(); j++)
				{
					images[i][j] = values[j]/norm;
				}
			}
		}


		void get_parallel(std::vector<std::vector<U>>& images, const std::vector<std::vector<U>>& str,const std::vector<std::vector<U>>& tangent, U k)
		{
			std::vector<std::vector<U>> amp(images.size(), std::vector<U>(2));
			for(size_t i = 1; i != images.size() - 1; i++)
			{
				U sum = 0.0;
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
				U diff = amp[i][1] - amp[i][0];
				for(size_t j = 0; j < images[i].size(); j++)
				{
					images[i][j] = k * diff * tangent[i][j];
				}
			}
		}
		U get_gradient(const vect& v, size_t n, U h = 0.000001) const
		{
			vect a = v;
			vect b = v;
			a[n]+=h;
			b[n]-=h;
			return (mainfunc(a) - mainfunc(b))/(2*h);
		}

		void get_gradients(std::vector<std::vector<U>>& images, const std::vector<std::vector<U>>& str)
		{
			for(size_t i = 1; i < images.size()-1; i++)
			{
				for(size_t j = 0; j < images[j].size(); j++)
				{
					images[i][j] = get_gradient(str[i],j,this->h);
				}
			}
		}
		void get_prepared(std::vector<std::vector<U>>& images, const std::vector<std::vector<U>>& tangent)
		{
			for(size_t i = 1; i < images.size()-1; i++)
			{
				U dot = scalar_product(images[i], tangent[i]);
				for(size_t j = 0; j < images[j].size(); j++)
				{
					images[i][j] = (images[i][j] - dot*tangent[i][j]);
				}
			}
		}
		void compute_mep()
		{
			path images(point_number+2, std::vector<U>(head.size()));

			for(size_t i=0; i<head.size(); i++)
			{
				images[0][i]=head[i];
				images[images.size()-1][i] = tail[i];
			}

			std::vector<std::vector<U>> tangents(point_number+2, std::vector<U>(head.size()));
			std::vector<std::vector<U>> springs(point_number+2, std::vector<U>(head.size()));
			std::vector<std::vector<U>> gradients(point_number+2, std::vector<U>(head.size()));
      
			std::vector<U> force(point_number*head.size(), 0);
			std::vector<U> old_force(point_number*head.size(), 0);
			std::vector<U> velocity(point_number*head.size(), 0);
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
        std::cout<<g<<std::endl;
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
        
        
        U sum = 0.0;
        U force_norm=norm(force);
        for(size_t i = 0; i < force.size(); i++)
						force[i] = force[i]/force_norm;  
				
        auto vf = scalar_product(velocity, force);
				if(vf > 0.0)
				{
					for(size_t i = 0; i < force.size(); i++)
						velocity[i] = vf*force[i];
				}
				else
				{
					for(size_t i = 0; i < velocity.size(); i++)
						velocity[i] = 0.0;
				}

				for(size_t i = 0; i < velocity.size(); i++)
				{
					velocity[i] += step*(old_force[i] + force[i]) / images[0].size();
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
					images[0] = findLocalMinimum(images[0],step,1);
					images[images.size()-1]=findLocalMinimum(images[images.size()-1],step,1);
				}
				else if(this->flex == 2)
				{
					images[images.size()-1]=findLocalMinimum(images[images.size()-1],step,1);

				}
				else if(this->flex == 3)
				{
					images[0] = findLocalMinimum(images[0],step,1);
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
  template <typename U>
	class STRING : public IMEP<U>
	{
    using path = std::vector<std::vector<U>>;
		using vect = std::vector<U>;
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
		U step;
		U h;
		U k;
		std::function<U(const vect&)> mainfunc;


	public:
		void set_parameters(const vect& _head,
							const vect& _tail,
							size_t pn = 10,
							U step = 0.01,
							size_t niter = 10000,
							U h = 1.0e-6,
							U k = 1.0)
		{

			this->head=_head;
			this->tail=_tail;
			this->point_number=pn;
			this->step=step;
			this->niter=niter;
			this->h=h;
			this->k=k;
		}
		void set_function(std::function<U(const vect&)> fn)
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
		std::vector<U> get_profile() const
		{
			return profile;
		}
		U get_rmse(const path& a, const path& b) const
		{
			if(a.size() != b.size())
			{
				/*std::cout<<"Вектора для подсчета RMSE разных размеров"<<std::endl;*/
				return -1;
			}
			else
			{
				U sum = 0.0;
				for(size_t i = 0; i < a.size(); i++)
				{
					for(size_t j = 0; j < a[i].size(); j++)
					{
						U diff = a[i][j] - b[i][j];
						sum += diff * diff;
					}
				}
				U mse = sum / (a.size()*a.front().size());
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

		void initialize_line(std::vector<std::vector<U>>& str, const std::vector<U>& head, const std::vector<U>& tail)
		{
			std::vector<U> es(head.size());
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
		std::vector<U> linear_interpolation(const std::vector<U>& gamma,
				const std::vector<U>& xp,
				const std::vector<U>& fp)
		{
			std::vector<U> interpolated_values;

			for(size_t i = 0; i < gamma.size(); i++)
			{
				auto lower = std::lower_bound(xp.begin(), xp.end(), gamma[i]);
				size_t idx = std::max(static_cast<int>(lower - xp.begin()) - 1, 0);

				U x1 = xp[idx];
				U x2 = xp[idx + 1];
				U y1 = fp[idx];
				U y2 = fp[idx + 1];

				U interpolated_y = y1 + (gamma[i] - x1) * (y2 - y1) / (x2 - x1);
				interpolated_values.push_back(interpolated_y);
			}

			return interpolated_values;
		}

		U norm(const std::vector<U>& v)
		{
			U sum = 0.0;
			for(U val : v)
			{
				sum += val * val;
			}
			return std::sqrt(sum);
		}


		U get_gradient(const vect& v, size_t n,U h = 0.0000001) const
		{
			vect a = v;
			vect b = v;
			a[n]+=h;
			b[n]-=h;
			return (mainfunc(a) - mainfunc(b))/(2*h);
		}

		std::vector<std::vector<U>> step_euler(std::vector<std::vector<U>>& string, U step)
		{

			path string_grad(string.size(),vect(string[0].size()));
			for(size_t i = 0; i < string.size(); i++)
			{
				for(size_t j=0; j< string[0].size(); j++)
					string_grad[i][j]=get_gradient(string[i],j);

			}
			U h = 0.0;
			for(size_t i = 0; i < string_grad.size(); i++)
			{
				U res=0;
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

		std::vector<U> linspace(U start, U end, int npts)
		{
			std::vector<U> result;
			U step = (end - start) / (npts - 1);
			for(int i = 0; i < npts; ++i)
			{
				result.push_back(start + i * step);
			}
			return result;
		}

		void compute_mep()
		{

			path string_xyz(point_number,std::vector<U>(head.size()));
			initialize_line(string_xyz,head,tail);
			std::vector<std::vector<U>> string;
			if(initial_path.empty()) { string=string_xyz; initial_path=string_xyz;}
			else   string = initial_path;
			mep_evolution.push_back(string);

			std::vector<std::vector<U>> old_string(string.size(), std::vector<U>(string[0].size(), 0.0));


			for(int i=1; i<niter+1; i++)
			{
        std::cout<<i<<std::endl;
				old_string=string;

				string=step_euler(string,step);

				std::vector<U> arclength;
				arclength.push_back(0.0);

				for(size_t it = 1; it < string.size(); it++)
				{
					std::vector<U> diff;
					for(size_t j = 0; j < string[it].size(); j++)
					{
						diff.push_back(string[it][j] - string[it - 1][j]);
					}
					U norm_diff = norm(diff);
					arclength.push_back(arclength.back() + norm_diff);
				}
				U arclength_last = arclength.back();
				for(U &value : arclength)
				{
					value /= arclength_last;
				}

				path xp_xyz(string[0].size(),std::vector<U>(arclength));
				std::vector<U> gamma=linspace(0,1,point_number);
				path fp_xyz(string[0].size(),std::vector<U>(arclength));
				path reparam_xyz(string[0].size(),std::vector<U>(string.size()));


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




#endif
