#ifndef LAZY_EXTENSIVE_UPDATER_HPP
#define LAZY_EXTENSIVE_UPDATER_HPP 1

#include "source/newtonian/two_dimensional/extensive_updater.hpp"

class LazyExtensiveUpdater: public ExtensiveUpdater
{
public:

  LazyExtensiveUpdater(const Tessellation& tess,
		       const PhysicalGeometry& pg);

  void operator()
  (const vector<Extensive>& fluxes,
   const PhysicalGeometry& pg,
   const Tessellation& tess,
   const double dt,
   const CacheData& cd,
   const vector<ComputationalCell>& cells,
   vector<Extensive>& extensives,
   double time) const;

private:
  const vector<double> lengths_;
};

#endif // LAZY_EXTENSIVE_UPDATER_HPP
