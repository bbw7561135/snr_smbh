#ifndef EDGE_LENGTH_CALCULATOR_HPP
#define EDGE_LENGTH_CALCULATOR_HPP 1

#include "source/misc/lazy_list.hpp"
#include "source/tessellation/tessellation.hpp"
#include "source/newtonian/two_dimensional/physical_geometry.hpp"

class EdgeLengthCalculator: public LazyList<double>
{
public:

  EdgeLengthCalculator(const Tessellation& tess,
		       const PhysicalGeometry& pg);

  size_t size(void) const;

  double operator[](size_t i) const;

private:
  const Tessellation& tess_;
  const PhysicalGeometry& pg_;
};

#endif // EDGE_LENGTH_CALCULATOR_HPP
