#include <vector>
#include <functional>
#include <cmath>

namespace mep
{
	struct IMEP
	{
		using path = std::vector<std::vector<double>>;
		using vect = std::vector<double>;
		virtual void set_parameters(const vect& _head, //neb
                                    const vect& _tail, 
                                    size_t pn = 10, 
                                    double step = 0.01, 
                                    size_t niter = 10000, 
                                    double h = 1.0e-6, 
                                    double k = 1.0) = 0;
		virtual void set_function(std::function<double(const vect&)> fn) = 0; //neb
		virtual double get_gradient(const vect& v, size_t n,double h) const = 0; 
		virtual vect get_profile() const =0;
		// 0 - свободные концы
		// 1 - фиксированные концы
		// 2 - фиксирована точка head, tail свободна
		// 3 - фиксирована точка tail, head свободна
        
		virtual void set_flexibility(size_t flex) = 0; 
		virtual void compute_mep() = 0; 
		virtual path get_mep() = 0; 
		virtual void set_optimal_path(const path& optimal_path) = 0; 
		virtual double get_rmse(const path& a, const path& b) const = 0; 
		virtual void set_initial_path(const path& a) = 0; 
    virtual path get_initial_path()const = 0;
		virtual vect get_rmse_evolution()  = 0;
		virtual std::vector<path> get_mep_evolution()  = 0; 
		
	};
}
