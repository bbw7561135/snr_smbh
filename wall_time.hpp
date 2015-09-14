#ifndef WALL_TIME_HPP
#define WALL_TIME_HPP 1

#include <ctime>
#include "source/newtonian/test_2d/main_loop_2d.hpp"

class WallTime: public DiagnosticFunction
{
public:

WallTime(const string& fname);

void operator()(const hdsim& sim);

private:
const string fname_;
const clock_t start_;
};

#endif // WALL_TIME_HPP
