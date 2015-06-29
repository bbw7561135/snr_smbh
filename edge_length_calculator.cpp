#include "edge_length_calculator.hpp"

EdgeLengthCalculator::EdgeLengthCalculator
(const Tessellation& tess,
 const PhysicalGeometry& pg):
  tess_(tess), pg_(pg) {}

size_t EdgeLengthCalculator::size(void) const
{
  return tess_.getAllEdges().size();
}

double EdgeLengthCalculator::operator[](size_t i) const
{
  return pg_.calcArea(tess_.getAllEdges()[i]);
}
