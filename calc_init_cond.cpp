#include "calc_init_cond.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/tessellation/shape_2d.hpp"

vector<ComputationalCell> calc_init_cond(const Tessellation& tess,
					 const Constants& c,
					 const EquationOfState& eos)
{
  vector<ComputationalCell> res(static_cast<size_t>(tess.GetPointNo()));
  for(size_t i=0;i<res.size();++i){
    const Vector2D r = tess.GetMeshPoint(static_cast<int>(i));
    res[i].density = 1.*c.proton_mass/pow(c.centi*c.meter,3);
    res[i].pressure = Uniform2D(1e-15)(r);
    res[i].velocity = r*c.wind_speed/abs(r);
    res[i].tracers["entropy"] = eos.dp2s
      (res[i].density, res[i].pressure);
    res[i].stickers["dummy"] = Circle(Vector2D(0,0),0.5*c.supernova_radius)(r);
  }
  return res;
}
