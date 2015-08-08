#include <ctime>
#include <fstream>
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include <fenv.h>
#include "constants.hpp"
#include "sim_data.hpp"
#include "my_main_loop.hpp"
#include "report_error.hpp"

using namespace std;

int main(void)
{
  const clock_t begin = clock();

  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

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

  const clock_t end = clock();
  ofstream f("wall_time.txt");
  f << static_cast<double>(end-begin)/CLOCKS_PER_SEC << endl;
  f.close();

  return 0;
}
