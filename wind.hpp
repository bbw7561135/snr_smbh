#ifndef WIND_HPP
#define WIND_HPP 1

#include "source/newtonian/two_dimensional/SourceTerm.hpp"

class Wind: public SourceTerm
{
public:

  Wind(const double& specific_mass_loss,
       const double& speed,
       const double& mass2thermal,
       const double& radius);

  vector<Extensive> operator()
  (const Tessellation& tess,
   const PhysicalGeometry& /*pg*/,
   const CacheData& cd,
   const vector<ComputationalCell>& /*cells*/,
   const vector<Extensive>& /*fluxes*/,
   const vector<Vector2D>& /*point_velocities*/,
   const double /*time*/) const;

private:
  const double specific_mass_loss_;
  const double speed_;
  const double mass2thermal_;
  const double radius_;
};

#endif // WIND_HPP
