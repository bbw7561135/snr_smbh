#ifndef SUPERNOA_HPP
#define SUPERNOA_HPP 1

#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/tessellation/shape_2d.hpp"

class Supernova: public Manipulate
{
public:

  Supernova(const Circle& hot_spot,
	    const double& mass,
	    const double& energy,
	    const double& time);

  void operator()(hdsim& sim);

private:
  const Circle hot_spot_;
  const double density_;
  const double pressure_;
  const double time_;
  mutable bool spent_;
};

#endif // SUPERNOA_HPP
