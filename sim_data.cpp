#include "sim_data.hpp"

SimData::SimData(const Constants& c):
  pg_(Vector2D(0,0), Vector2D(0,1)),
  outer_(c.lower_left, c.upper_right),
  init_points_(clip_grid
	       (RightRectangle(c.lower_left+Vector2D(0.001,0), c.upper_right),
		complete_grid(0.1*c.parsec,
			      abs(c.upper_right-c.lower_left),
			      0.001))),
  tess_(init_points_, outer_),
  eos_(c.adiabatic_index),
  rs_(),
  //    raw_point_motion_(),
  //    point_motion_(raw_point_motion_,eos_),
  alt_point_motion_(),
  gravity_acc_(c.gravitation_constant*c.black_hole_mass,
	       0.001*c.parsec,
	       c.parsec*Vector2D(0,0)),
  gravity_force_(gravity_acc_),
  geom_force_(pg_.getAxis()),
  wind_(1e-3*c.solar_mass/c.year/(4.*M_PI*pow(0.4*c.parsec,3)/3.),
	c.wind_speed,
	c.boltzmann_constant*1e4/(5./3.-1)/c.proton_mass,
	0.4*c.parsec),
  force_(VectorInitialiser<SourceTerm*>(&gravity_force_)(&wind_)(&geom_force_)()),
  //    force_(VectorInitializer<SourceTerm*>(&gravity_force_)()),
  tsf_(0.3, "dt_log.txt",c.gravitation_constant*c.black_hole_mass),
  fc_(rs_),
  eu_(tess_,pg_),
  cu_(),
  sim_(tess_,
       outer_,
       pg_,
       calc_init_cond(tess_,c,eos_),
       eos_,
       alt_point_motion_,
       force_,
       tsf_,
       fc_,
       eu_,
       cu_) {}

hdsim& SimData::getSim(void)
{
  return sim_;
}
