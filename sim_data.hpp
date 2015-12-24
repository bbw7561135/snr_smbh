#ifndef SIM_DATA_HPP
#define SIM_DATA_HPP 1

#include "constants.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
//#include "source/tessellation/static_voronoi_mesh.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/newtonian/two_dimensional/source_terms/CenterGravity.hpp"
#include "source/newtonian/two_dimensional/source_terms/cylindrical_complementary.hpp"
#include "wind.hpp"
#include "source/newtonian/two_dimensional/source_terms/SeveralSources.hpp"
#include "log_cfl.hpp"
#include "sink_flux.hpp"
#include "lazy_extensive_updater.hpp"
#include "entropy_fix.hpp"
#include "source/tessellation/right_rectangle.hpp"
#include "complete_grid.hpp"
#include "source/newtonian/test_2d/clip_grid.hpp"
#include "source/misc/vector_initialiser.hpp"
#include "calc_init_cond.hpp"
#include "source/newtonian/two_dimensional/stationary_box.hpp"

class SimData
{

public:

  SimData(const Constants& c);

  hdsim& getSim(void);
  
private:
  const CylindricalSymmetry pg_;
  const SquareBox outer_;
  const vector<Vector2D> init_points_;
  //StaticVoronoiMesh tess_;
  VoronoiMesh tess_;
  const IdealGas eos_;
  const Hllc rs_;
  //  Lagrangian raw_point_motion_;
  //  RoundCells point_motion_;
  Eulerian alt_point_motion_;
  const StationaryBox evc_;
  CenterGravity gravity_acc_;
  ConservativeForce gravity_force_;
  CylindricalComplementary geom_force_;
  Wind wind_;
  //  RadiativeCooling rad_cool_;
  SeveralSources force_;
  //  const SimpleCFL tsf_;
  const LogCFL tsf_;
  const SinkFlux fc_;
  //const SimpleCellUpdater cu_;
  const LazyExtensiveUpdater eu_;
  const EntropyFix cu_;
  hdsim sim_;
};

#endif // SIM_DATA_HPP
