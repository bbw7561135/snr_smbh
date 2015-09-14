#include "wall_time.hpp"
#include "source/misc/simple_io.hpp"

WallTime::WallTime(const string& fname):
  fname_(fname), start_(clock()) {}

void WallTime::operator()(const hdsim& /*sim*/)
{
  write_number
    (static_cast<double>(clock()-start_)/CLOCKS_PER_SEC,
     fname_);
}
