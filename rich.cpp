#include <iostream>
#include <cmath>
#include <limits>
#include "source/tessellation/static_voronoi_mesh.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/newtonian/two_dimensional/interpolations/pcm2d.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/newtonian/two_dimensional/source_terms/cylindrical_complementary.hpp"
#include "source/newtonian/two_dimensional/source_terms/CenterGravity.hpp"
#include "source/newtonian/two_dimensional/source_terms/SeveralSources.hpp"
#include "source/misc/simple_io.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/newtonian/test_2d/kill_switch.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/tessellation/RoundGrid.hpp"
#include "source/misc/int2str.hpp"
#include "source/newtonian/two_dimensional/interpolations/linear_gauss_consistent.hpp"
#include "source/newtonian/test_2d/consecutive_snapshots.hpp"
#include <boost/foreach.hpp>
#include "source/misc/utils.hpp"
#include "source/tessellation/shape_2d.hpp"
#include "source/newtonian/test_2d/piecewise.hpp"
#include "source/mpi/MeshPointsMPI.hpp"
#include "source/mpi/mpi_macro.hpp"
#include "source/mpi/ConstNumberPerProc.hpp"
#include "source/misc/hdf5_utils.hpp"
#include "source/newtonian/test_2d/contour.hpp"
#include "source/newtonian/two_dimensional/custom_evolutions/Ratchet.cpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/newtonian/two_dimensional/simple_flux_calculator.hpp"
#include "source/newtonian/two_dimensional/simple_cell_updater.hpp"
#include "source/misc/vector_initialiser.hpp"
#include "source/misc/horner.hpp"
#include <fenv.h>
#include "write_conserved.hpp"
#include "constants.hpp"
#include "wind.hpp"
#include "supernova.hpp"
#include "source/tessellation/right_rectangle.hpp"
#include "source/newtonian/test_2d/clip_grid.hpp"
#include "complete_grid.hpp"
#include "sink_flux.hpp"
#include "calc_init_cond.hpp"
#include "log_cfl.hpp"
#include "entropy_fix.hpp"
#include "edge_length_calculator.hpp"

using namespace std;
using namespace simulation2d;

namespace {

  bool bracketed(double low, double arg, double high)
  {
    return arg>=low and high>arg;
  }

  class LazyExtensiveUpdater: public ExtensiveUpdater
  {
  public:

    LazyExtensiveUpdater(const Tessellation& tess,
			 const PhysicalGeometry& pg):
      lengths_(serial_generate(EdgeLengthCalculator(tess,pg))) {}

    void operator()
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

  private:
    const vector<double> lengths_;
  };
}

class SimData
{

public:

  SimData(const Constants& c):
   pg_(Vector2D(0,0), Vector2D(0,1)),
    outer_(c.lower_left, c.upper_right),
    init_points_(clip_grid
		 (RightRectangle(c.lower_left+Vector2D(0.001,0), c.upper_right),
		  complete_grid(0.1*c.parsec,
				abs(c.upper_right-c.lower_left),
				0.005*2))),
    tess_(init_points_, outer_),
    eos_(c.adiabatic_index),
    rs_(),
    //    raw_point_motion_(),
    //    point_motion_(raw_point_motion_,eos_),
    alt_point_motion_(),
    gravity_acc_(c.gravitation_constant*c.black_hole_mass,
		 0.001*c.parsec,
		 c.parsec*Vector2D(0,0)),
    gravity_force_(gravity_acc_),
    geom_force_(pg_.getAxis()),
    wind_(1e-3*c.solar_mass/c.year/(4.*M_PI*pow(0.4*c.parsec,3)/3.),
	  1e3*c.kilo*c.meter/c.second,
	  c.boltzmann_constant*1e4/(5./3.-1)/c.proton_mass,
	  0.4*c.parsec),
    force_(VectorInitialiser<SourceTerm*>(&gravity_force_)(&wind_)(&geom_force_)()),
  //    force_(VectorInitializer<SourceTerm*>(&gravity_force_)()),
    tsf_(0.3, "dt_log.txt",c.gravitation_constant*c.black_hole_mass),
    fc_(rs_),
    eu_(tess_,pg_),
    cu_(),
    sim_(tess_,
	 outer_,
	 pg_,
	 calc_init_cond(tess_,c,eos_),
	 eos_,
	 alt_point_motion_,
	 force_,
	 tsf_,
	 fc_,
	 eu_,
	 cu_) {}

  hdsim& getSim(void)
  {
    return sim_;
  }
  
private:
  const CylindricalSymmetry pg_;
  const SquareBox outer_;
  const vector<Vector2D> init_points_;
  StaticVoronoiMesh tess_;
  const IdealGas eos_;
  const Hllc rs_;
  //  Lagrangian raw_point_motion_;
  //  RoundCells point_motion_;
  Eulerian alt_point_motion_;
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

class CustomDiagnostic
{
public:

  CustomDiagnostic(double dt,
		   double t_start):
    dt_(dt),
    t_next_(t_start),
    counter_(0) {}

  void operator()(const hdsim& sim)
  {
    if(sim.getTime()>t_next_){
      write_snapshot_to_hdf5(sim, "snapshot_"+int2str(counter_)+".h5");

      t_next_ += dt_;
      ++counter_;
    }
  }

private:
  const double dt_;
  double t_next_;
  int counter_;
};

namespace {
  class MultipleDiagnostics: public DiagnosticFunction
  {
  public:

    MultipleDiagnostics(const vector<DiagnosticFunction*>& diag_list):
      diag_list_(diag_list) {}

    void operator()(const hdsim& sim)
    {
      BOOST_FOREACH(DiagnosticFunction* df, diag_list_)
	{
	  (*df)(sim);
	}
    }

  private:
    const vector<DiagnosticFunction*> diag_list_;
  };
}

namespace {

  class WriteCycle: public DiagnosticFunction
  {
  public:

    WriteCycle(const string& fname):
      fname_(fname) {}

    void operator()(const hdsim& sim)
    {
      write_number(sim.getCycle(),fname_);
    }

  private:
    const string fname_;
  };

  void my_main_loop(hdsim& sim, const Constants& c)
  {
    const double tf = 1e3*c.year;
    SafeTimeTermination term_cond_raw(tf,1e6);
    ConsecutiveSnapshots diag1(new ConstantTimeInterval(tf/10),
			       new Rubric("snapshot_",".h5"));
    WriteTime diag2("time.txt");
    WriteConserved diag4("conserved.txt");
    //Bremsstrahlung diag5("luminosity_history.txt",c);
    WriteCycle diag6("cycle.txt");
    vector<DiagnosticFunction*> diag_list;
    diag_list.push_back(&diag1);
    diag_list.push_back(&diag2);
    //    diag_list.push_back(&diag3);
    diag_list.push_back(&diag4);
    //    diag_list.push_back(&diag5);
    diag_list.push_back(&diag6);
    MultipleDiagnostics diag(diag_list);
    const Circle hot_spot(Vector2D(0,-c.offset),
			  c.supernova_radius);
    Supernova manip(hot_spot,
		    c.supernova_mass,
		    c.supernova_energy,
		    2e4*c.year);
    main_loop(sim, 
	      term_cond_raw,
	      &hdsim::TimeAdvance, 
	      &diag,
	      &manip);
  }
}

namespace {
  void report_error(UniversalError const& eo)
  {
    cout << eo.GetErrorMessage() << endl;
    cout.precision(14);
    for(size_t i=0;i<eo.GetFields().size();++i)
      cout << eo.GetFields()[i] << " = "
	   << eo.GetValues()[i] << endl;
  }
}

int main(void)
{
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  //  MPI_Init(NULL, NULL);

  assert(abs(horner(VectorInitialiser<double>(1)(2)(3)(),2)-11)<1e-6);

  try{

    Constants c;

    SimData sim_data(c);
    hdsim& sim = sim_data.getSim();

    write_snapshot_to_hdf5(sim,
			   "initial.h5");

    my_main_loop(sim,c);

    write_snapshot_to_hdf5(sim, 
			   "final.h5");

  }
  catch(const UniversalError& eo){
    report_error(eo);
    throw;
  }

  return 0;
}
