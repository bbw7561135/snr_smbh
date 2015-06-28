#include <cmath>
#include "complete_grid.hpp"
#include "logarithmic_spiral.hpp"
#include "hexagonal_grid.hpp"
#include "source/misc/utils.hpp"

vector<Vector2D> complete_grid(double r_inner,
			       double r_outer,
			       double alpha)
{
  const vector<Vector2D> inner = 
    centered_hexagonal_grid(r_inner*alpha*2*M_PI,
			    r_inner);
  const vector<Vector2D> outer =
    centered_logarithmic_spiral(r_inner,
				r_outer,
				alpha,
				Vector2D(0,0));
  return join(inner, outer);
}
