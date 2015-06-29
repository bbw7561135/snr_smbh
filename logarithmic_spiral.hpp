#ifndef LOGARITHMIC_SPIRAL_HPP
#define LOGARITHMIC_SPIRAL_HPP 1

#include <vector>
#include "source/tessellation/geometry.hpp"

using std::vector;

vector<Vector2D> centered_logarithmic_spiral(double r_min,
					     double r_max,
					     double alpha,
					     const Vector2D& center);

#endif // LOGARITHMIC_SPIRAL_HPP
