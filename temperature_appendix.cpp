#include "temperature_appendix.hpp"

TemperatureAppendix::TemperatureAppendix(double m2k): m2k_(m2k) {}

string TemperatureAppendix::getName(void) const
{
  return "temperature";
}

vector<double> TemperatureAppendix::operator()(const hdsim& sim) const
{
  const vector<ComputationalCell>& cells = sim.getAllCells();
  vector<double> res(cells.size());
  for(size_t i=0;i<res.size();++i){
    const ComputationalCell& cell = cells[i];
    res[i] = m2k_*cell.pressure/cell.density;
  }
  return res;
}
