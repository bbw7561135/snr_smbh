#include <fstream>
#include "write_conserved.hpp"

using namespace std;

WriteConserved::WriteConserved(const string& fname):
  time_cons_(), fname_(fname) {}

void WriteConserved::operator()(const hdsim& sim)
{
  Extensive buf;
  buf.mass = 0;
  buf.momentum.x = 0;
  buf.momentum.y = 0;
  buf.energy = 0;

  for(size_t i=0;i<sim.getAllExtensives().size();++i){
    if(sim.getAllCells()[i].stickers.find("dummy")->second){
      const Extensive& temp = sim.getAllExtensives()[i];
      buf.mass += temp.mass;
      buf.momentum += temp.momentum;
      buf.energy += temp.energy;
    }
  }
  time_cons_.push_back(pair<double,Extensive>
		       (sim.getTime(),buf));
}

WriteConserved::~WriteConserved(void)
{
  ofstream f(fname_.c_str());
  for(size_t i=0;i<time_cons_.size();++i)
    f << time_cons_[i].first << " "
      << time_cons_[i].second.mass << " "
      << time_cons_[i].second.momentum.x << " "
      << time_cons_[i].second.momentum.y << " "
      << time_cons_[i].second.energy << endl;
  f.close();
}
