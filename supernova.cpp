#include "supernova.hpp"

namespace {
  double calc_sphere_volume(double radius)
  {
    return 4.0*M_PI*pow(radius,3)/3.;
  }
}

Supernova::Supernova(const Circle& hot_spot,
		     const double& mass,
		     const double& energy,
		     const double& time):
  hot_spot_(hot_spot),
  density_(mass/calc_sphere_volume(hot_spot.getRadius())),
  pressure_(energy/calc_sphere_volume(hot_spot.getRadius())*(5./3.-1)),
  time_(time),
  spent_(false) {}

void Supernova::operator()(hdsim& sim)
{
  if(time_>sim.getTime() || spent_)
    return;
  spent_ = true;
  vector<ComputationalCell>& cells = sim.getAllCells();
  for(size_t i=0;i<cells.size();++i){
    const Vector2D r = sim.getTessellation().GetMeshPoint(static_cast<int>(i));
    if(hot_spot_(r)){
      ComputationalCell& cell = cells[i]; 
      const double old_density = cell.density;
      cell.density += density_;
      cell.pressure += pressure_;
      cell.velocity = cell.velocity*old_density/cell.density;
    }
  }
  sim.recalculateExtensives();
}
