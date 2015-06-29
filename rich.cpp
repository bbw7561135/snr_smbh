#include <iostream>
#include <cmath>
#include <limits>
#include "source/newtonian/two_dimensional/interpolations/pcm2d.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/hydro_boundary_conditions/RigidWallHydro.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
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
#include "source/misc/horner.hpp"
#include <fenv.h>
#include "write_conserved.hpp"
#include "constants.hpp"
#include "supernova.hpp"
#include "edge_length_calculator.hpp"
#include "sim_data.hpp"
#include "source/newtonian/test_2d/multiple_diagnostics.hpp"
#include "write_cycle.hpp"

using namespace std;
using namespace simulation2d;

namespace {

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
