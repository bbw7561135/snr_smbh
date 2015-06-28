#include "supernova.hpp"

Supernova::Supernova(const Circle& hot_spot,
		     const double& mass,
		     const double& energy,
		     const double& time):
  hot_spot_(hot_spot),
  density_(mass/pow(hot_spot.getRadius(),3)),
  pressure_(energy/pow(hot_spot.getRadius(),3)),
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
      cell.density = density_;
      cell.pressure = pressure_;
    }
  }
  sim.recalculateExtensives();
}
