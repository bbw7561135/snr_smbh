#include "write_cycle.hpp"

WriteCycle::WriteCycle(const string& fname):
  fname_(fname) {}

void WriteCycle::operator()(const hdsim& sim)
{
  write_number(sim.getCycle(),fname_);
}
