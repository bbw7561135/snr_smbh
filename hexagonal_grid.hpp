#ifndef HEXAGONAL_GRID_HPP
#define HEXAGONAL_GRID_HPP 1

#include <vector>
#include "source/tessellation/geometry.hpp"

using std::vector;

vector<Vector2D> centered_hexagonal_grid
(double r_min, double r_max);

#endif // HEXAGONAL_GRID_HPP
