#include <cmath>
#include "logarithmic_spiral.hpp"
#include "source/misc/utils.hpp"

vector<Vector2D> centered_logarithmic_spiral(double r_min,
					     double r_max,
					     double alpha,
					     const Vector2D& center)
{
  const double theta_max = log(r_max/r_min)/alpha;
  const vector<double> theta_list = 
    arange(0,theta_max,2*M_PI*alpha/(1-0.5*alpha));
  
  vector<double> r_list(theta_list.size(),0);
  for(size_t i=0;i<r_list.size();++i)
    r_list.at(i) = r_min*exp(alpha*theta_list.at(i));
  
  vector<Vector2D> res(r_list.size());
  for(size_t i=0;i<res.size();++i)
    res[i] = center+r_list[i]*Vector2D(cos(theta_list.at(i)),
				       sin(theta_list.at(i)));
  return res;
}
