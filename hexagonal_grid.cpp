#include <cmath>
#include "hexagonal_grid.hpp"
#include "source/misc/utils.hpp"

using std::size_t;
using std::max;

vector<Vector2D> centered_hexagonal_grid(double r_min,
					 double r_max)
{
  const vector<double> r_list = arange(0,r_max,r_min);
  vector<Vector2D> res;
  for(size_t i=0;i<r_list.size();++i){
    const size_t angle_num = max<size_t>(6*i,1);
    vector<double> angle_list(angle_num,0);
    for(size_t j=0;j<angle_num;++j)
      angle_list.at(j) = 2*M_PI*static_cast<double>(j)/static_cast<double>(angle_num);
    for(size_t j=0;j<angle_num;++j)
      res.push_back(r_list.at(i)*Vector2D(cos(angle_list.at(j)),
					  sin(angle_list.at(j))));
  }
  return res;
}
