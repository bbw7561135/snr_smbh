#ifndef LOG_CFL_HPP
#define LOG_CFL_HPP 1

#include "source/newtonian/two_dimensional/time_step_function.hpp"

class LogCFL: public TimeStepFunction
{
public:

  LogCFL(double cfl, const string& fname, double gm);

  double operator()
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
   const EquationOfState& eos,
   const vector<Vector2D>& point_velocities,
   const double /*time*/) const;

private:
  const double cfl_;
  const string fname_;
  const double gm_;
};

#endif // LOG_CFL_HPP
