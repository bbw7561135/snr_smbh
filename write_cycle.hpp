#ifndef WRITE_CYCLE_HPP
#define WRITE_CYCLE_HPP 1

#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/misc/simple_io.hpp"

class WriteCycle: public DiagnosticFunction
{
public:

  WriteCycle(const string& fname);

  void operator()(const hdsim& sim);

private:
  const string fname_;
};

#endif // WRITE_CYCLE_HPP
