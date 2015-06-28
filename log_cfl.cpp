#include <cmath>
#include <fstream>
#include "log_cfl.hpp"
#include "source/misc/lazy_list.hpp"

using std::ofstream;
using std::endl;

namespace {
  class TimeStepCalculator: public LazyList<double>
  {
  public:

    TimeStepCalculator(const Tessellation& tess,
		       const vector<ComputationalCell>& cells,
		       const EquationOfState& eos,
		       const vector<Vector2D>& point_velocities,
		       const double gm):
      tess_(tess), cells_(cells), 
      point_velocities_(point_velocities), eos_(eos), gm_(gm) {}

    size_t size(void) const 
    {
      return cells_.size();
    }

    double operator[](size_t i) const
    {
      const double kepler_velocity = sqrt(gm_/abs(tess_.GetMeshPoint(static_cast<int>(i))));
      return tess_.GetWidth(static_cast<int>(i))/
	(eos_.dp2c(cells_[i].density, cells_[i].pressure)+
	 abs(cells_[i].velocity)+
	 kepler_velocity+
	 abs(point_velocities_[i]));
    }

  private:
    const Tessellation& tess_;
    const vector<ComputationalCell>& cells_;
    const vector<Vector2D>& point_velocities_;
    const EquationOfState& eos_;
    const double gm_;
  };
}
  

LogCFL::LogCFL(double cfl, const string& fname, double gm):
  cfl_(cfl), fname_(fname), gm_(gm) {}

double LogCFL::operator()
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
   const EquationOfState& eos,
   const vector<Vector2D>& point_velocities,
   const double /*time*/) const
{
  const TimeStepCalculator tsc(tess,cells,eos,point_velocities,gm_);
  double res = 10000;
  size_t argmin = 0;
  for(size_t i=0;i<tsc.size();++i){
    const double candidate = tsc[i];
    if(candidate<res && !cells[i].stickers.find("dummy")->second){
      res = candidate;
      argmin = i;
    }	
  }
  ofstream f(fname_.c_str());
  f << argmin << endl;
  f << tess.GetMeshPoint(static_cast<int>(argmin)).x << ", ";
  f << tess.GetMeshPoint(static_cast<int>(argmin)).y << endl;
  f << tess.GetWidth(static_cast<int>(argmin)) << endl;
  f << eos.dp2c(cells[argmin].density, cells[argmin].pressure) << endl;
  f << cells[argmin].velocity.x << ", ";
  f << cells[argmin].velocity.y << endl;
  f << point_velocities[argmin].x << ", ";
  f << point_velocities[argmin].y << endl;
  f.close();
  return cfl_*res;
}
