#ifndef COMPLETE_GRIDH_HPP
#define COMPLETE_GRIDH_HPP 1

#include <vector>
#include "source/tessellation/geometry.hpp"

using std::vector;

vector<Vector2D> complete_grid(double r_inner,
			       double r_outer,
			       double alpha);

#endif // COMPLETE_GRIDH_HPP
