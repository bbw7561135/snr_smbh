#ifndef SINK_FLUX_HPP
#define SINK_FLUX_HPP 1

#include "source/newtonian/two_dimensional/flux_calculator_2d.hpp"
#include "source/newtonian/common/riemann_solver.hpp"
#include "source/newtonian/two_dimensional/extensive.hpp"
#include "source/tessellation/tessellation.hpp"

using std::vector;

class SinkFlux: public FluxCalculator
{
public:

  SinkFlux(const RiemannSolver& rs);

  vector<Extensive> operator()
  (const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const vector<Extensive>& extensives,
   const CacheData& cd,
   const EquationOfState& eos,
   const double /*time*/,
   const double /*dt*/) const;

private:
  const RiemannSolver& rs_;

  const Conserved calcHydroFlux
  (const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const EquationOfState& eos,
   const size_t i) const;
};

#endif // SINK_FLUX_HPP
