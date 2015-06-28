#ifndef ENTROPY_FIX_HPP
#define ENTROPY_FIX_HPP 1

#include "source/newtonian/two_dimensional/cell_updater_2d.hpp"

class EntropyFix: public CellUpdater
{
public:
  EntropyFix(const double thres=1e-9);

  vector<ComputationalCell> operator()
  (const Tessellation& /*tess*/,
   const PhysicalGeometry& /*pg*/,
   const EquationOfState& eos,
   const vector<Extensive>& extensives,
   const vector<ComputationalCell>& old,
   const CacheData& cd) const;

private:
  const double thres_;
};

#endif // ENTROPY_FIX_HPP
