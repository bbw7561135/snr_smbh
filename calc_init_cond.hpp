#ifndef CALC_INIT_COND_HPP
#define CALC_INIT_COND_HPP 1

#include <vector>
#include <cmath>
#include "source/newtonian/two_dimensional/computational_cell_2d.hpp"
#include "source/tessellation/tessellation.hpp"
#include "constants.hpp"
#include "source/newtonian/common/equation_of_state.hpp"

using std::vector;

vector<ComputationalCell> calc_init_cond(const Tessellation& tess,
					 const Constants& c,
					 const EquationOfState& eos);

#endif // CALC_INIT_COND_HPP
