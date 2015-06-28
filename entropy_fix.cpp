#include "entropy_fix.hpp"

EntropyFix::EntropyFix(const double thres):
  thres_(thres) {}

vector<ComputationalCell> EntropyFix::operator()
  (const Tessellation& /*tess*/,
   const PhysicalGeometry& /*pg*/,
   const EquationOfState& eos,
   const vector<Extensive>& extensives,
   const vector<ComputationalCell>& old,
   const CacheData& cd) const
{
  vector<ComputationalCell> res = old;
  for(size_t i=0;i<extensives.size();++i){
    assert(extensives[i].tracers.find("entropy")->second>0);
    if(old[i].stickers.find("dummy")->second)
      continue;
    /*
      const double volume = pg.calcVolume
      (serial_generate(CellEdgesGetter(tess,static_cast<int>(i))));
    */
    const double volume = cd.volumes[i];
    res[i].density = extensives[i].mass/volume;
    assert(res[i].density>0);
    res[i].velocity = extensives[i].momentum / extensives[i].mass;
    const double total_energy = extensives[i].energy/extensives[i].mass;
    const double kinetic_energy = 
      0.5*ScalarProd(res[i].velocity, res[i].velocity);
    const double energy = 
      total_energy - kinetic_energy;
    for(std::map<std::string,double>::const_iterator it =
	  extensives[i].tracers.begin();
	it!=extensives[i].tracers.end();++it)
      res[i].tracers[it->first] = it->second/extensives[i].mass;
    //	assert(energy>thres_*kinetic_energy);
    if(energy>thres_*kinetic_energy){
      res[i].pressure = eos.de2p(res[i].density, energy);
      res[i].tracers["entropy"] = eos.dp2s
	(res[i].density, res[i].pressure);
      assert(res[i].pressure>0);
      //	  assert(res[i].tracers["entropy"]/old[i].tracers.find("entropy")->second<1000);
    }
    else{
      assert(res[i].tracers.count("entropy")>0);
      res[i].pressure = eos.sd2p(res[i].tracers.find("entropy")->second,
				 res[i].density);
      assert(res[i].pressure>0);
    }
  }
  return res;
}
