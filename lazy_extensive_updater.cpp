#include "lazy_extensive_updater.hpp"
#include "edge_length_calculator.hpp"

namespace {
  bool bracketed(double low, double arg, double high)
  {
    return arg>=low and high>arg;
  }
}

LazyExtensiveUpdater::LazyExtensiveUpdater
(const Tessellation& tess,
 const PhysicalGeometry& pg):
  lengths_(serial_generate(EdgeLengthCalculator(tess,pg))) {}

void LazyExtensiveUpdater::operator()
  (const vector<Extensive>& fluxes,
   const PhysicalGeometry& /*pg*/,
   const Tessellation& tess,
   const double dt,
   const CacheData& cd,
   const vector<ComputationalCell>& cells,
   vector<Extensive>& extensives) const
{
  const vector<Edge>& edge_list = tess.getAllEdges();
  for(size_t i=0;i<cells.size();++i){
    extensives[i].tracers["entropy"] = 
      cd.volumes[i]*cells[i].density*
      cells[i].tracers.find("entropy")->second;
  }
  for(size_t i=0;i<edge_list.size();++i){
    const Edge& edge = edge_list[i];
    const Extensive delta = dt*lengths_[i]*fluxes[i];
    if(bracketed(0,edge.neighbors.first,tess.GetPointNo())){
      extensives[static_cast<size_t>(edge.neighbors.first)] -=
	delta;
      assert(extensives[static_cast<size_t>(edge.neighbors.first)].tracers.find("entropy")->second>0);
    }
    if(bracketed(0,edge.neighbors.second,tess.GetPointNo())){
      extensives[static_cast<size_t>(edge.neighbors.second)] +=
	delta;
      assert(extensives[static_cast<size_t>(edge.neighbors.second)].tracers.find("entropy")->second>0);
    }
  }
}
