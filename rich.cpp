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
//#include "source/newtonian/two_dimensional/hydro_boundary_conditions/FreeFlow.hpp"
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

using namespace std;
using namespace simulation2d;

namespace {

  vector<Vector2D> centered_hexagonal_grid(double r_min,
					   double r_max)
  {
    const vector<double> r_list = arange(0,r_max,r_min);
    vector<Vector2D> res;
    for(size_t i=0;i<r_list.size();++i){
      const size_t angle_num = max<size_t>(6*i,1);
      vector<double> angle_list(angle_num,0);
      for(size_t j=0;j<angle_num;++j)
	angle_list.at(j) = 2*M_PI*static_cast<double>(j)/static_cast<double>(angle_num);
      for(size_t j=0;j<angle_num;++j)
	res.push_back(r_list.at(i)*Vector2D(cos(angle_list.at(j)),
					    sin(angle_list.at(j))));
    }
    return res;
  }

  vector<Vector2D> centered_logarithmic_spiral(double r_min,
					       double r_max,
					       double alpha,
					       const Vector2D& center)
  {
    const double theta_max = log(r_max/r_min)/alpha;
    const vector<double> theta_list = 
      arange(0,theta_max,2*M_PI*alpha/(1-0.5*alpha));
  
    vector<double> r_list(theta_list.size(),0);
    for(size_t i=0;i<r_list.size();++i)
      r_list.at(i) = r_min*exp(alpha*theta_list.at(i));
  
    vector<Vector2D> res(r_list.size());
    for(size_t i=0;i<res.size();++i)
      res[i] = center+r_list[i]*Vector2D(cos(theta_list.at(i)),
					 sin(theta_list.at(i)));
    return res;
  }

  vector<Vector2D> complete_grid(double r_inner,
				 double r_outer,
				 double alpha)
  {
    const vector<Vector2D> inner = 
      centered_hexagonal_grid(r_inner*alpha*2*M_PI,
			      r_inner);
    const vector<Vector2D> outer =
      centered_logarithmic_spiral(r_inner,
				  r_outer,
				  alpha,
				  Vector2D(0,0));
    return join(inner, outer);
  }

  /*
  template<class T> class VectorInitializer
  {
  public:

    VectorInitializer(const T& t):
      buf_(1,t) {}

    VectorInitializer& operator()(const T& t)
    {
      buf_.push_back(t);
      return *this;
    }

    vector<T> operator()(void)
    {
      return buf_;
    }

  private:
    vector<T> buf_;
  };
  */

  vector<ComputationalCell> calc_init_cond(const Tessellation& tess,
					   const Constants& c,
					   const EquationOfState& eos)
  {
    vector<ComputationalCell> res(static_cast<size_t>(tess.GetPointNo()));
    for(size_t i=0;i<res.size();++i){
      const Vector2D r = tess.GetMeshPoint(static_cast<int>(i));
      res[i].density = 1.*c.proton_mass/pow(c.centi*c.meter,3);
      res[i].pressure = Uniform2D(1e-15)(r);
      res[i].velocity = Vector2D(0,0);
      //      res[i].velocity = 1000*c.kilo*c.meter/c.second*r/abs(r);
      /*
      res[i].tracers["ejecta"] = Piecewise
	(Circle(Vector2D(0,-c.offset),c.supernova_radius),
	 Uniform2D(1), Uniform2D(0))(r);
      */
      res[i].tracers["entropy"] = eos.dp2s
	(res[i].density, res[i].pressure);
      res[i].stickers["dummy"] = Circle(Vector2D(0,0),0.5*c.supernova_radius)(r);
    }
    return res;
  }

  double calc_tracer_flux(const Edge& edge,
			  const Tessellation& tess,
			  const vector<ComputationalCell>& cells,
			  const string& name,
			  const Conserved& hf)
  {
    if(hf.Mass>0 && 
       edge.neighbors.first>0 && 
       edge.neighbors.first<tess.GetPointNo())
      return hf.Mass*
	cells[static_cast<size_t>(edge.neighbors.first)].tracers.find(name)->second;
    if(hf.Mass<0 && 
       edge.neighbors.second>0 && 
       edge.neighbors.second<tess.GetPointNo())
      return hf.Mass*
	cells[static_cast<size_t>(edge.neighbors.second)].tracers.find(name)->second;
    return 0;    
  }

  class SinkFlux: public FluxCalculator
  {
  public:

    SinkFlux(const RiemannSolver& rs):
      rs_(rs) {}

    vector<Extensive> operator()
    (const Tessellation& tess,
     const vector<Vector2D>& point_velocities,
     const vector<ComputationalCell>& cells,
     const EquationOfState& eos,
     const double /*time*/,
     const double /*dt*/) const
    {
      vector<Extensive> res(tess.getAllEdges().size());
      for(size_t i=0;i<tess.getAllEdges().size();++i){
	const Conserved hydro_flux =
	  calcHydroFlux(tess, point_velocities,
			cells, eos, i);
	res[i].mass = hydro_flux.Mass;
	res[i].momentum = hydro_flux.Momentum;
	res[i].energy = hydro_flux.Energy;
	for(std::map<std::string,double>::const_iterator it =
	      cells.front().tracers.begin();
	    it!=cells.front().tracers.end();++it)
	  res[i].tracers[it->first] =
	    calc_tracer_flux(tess.getAllEdges()[i],
			     tess,cells,it->first,hydro_flux);
	     
      }
      return res;
    }

  private:
    const RiemannSolver& rs_;

    const Conserved calcHydroFlux
    (const Tessellation& tess,
     const vector<Vector2D>& point_velocities,
     const vector<ComputationalCell>& cells,
     const EquationOfState& eos,
     const size_t i) const
    {
      const Edge& edge = tess.GetEdge(static_cast<int>(i));
      const std::pair<bool,bool> flags
	(edge.neighbors.first>=0 && edge.neighbors.first<tess.GetPointNo(),
	 edge.neighbors.second>=0 && edge.neighbors.second<tess.GetPointNo());
      assert(flags.first || flags.second);
      if(!flags.first){
	const size_t right_index = 
	  static_cast<size_t>(edge.neighbors.second);
	const ComputationalCell& right_cell = cells[right_index];
	if(right_cell.stickers.find("dummy")->second)
	  return Conserved();
	const Vector2D p = Parallel(edge);
	const Primitive right = convert_to_primitive(right_cell,eos);
	//	const Primitive left = reflect(right,p);
	const Primitive left = right;
	const Vector2D n = remove_parallel_component
	  (tess.GetMeshPoint(edge.neighbors.second) - 
	   edge.vertices.first, p);
	return rotate_solve_rotate_back
	  (rs_, left, right, 0, n, p);
      }
      if(!flags.second){
	const size_t left_index = 
	  static_cast<size_t>(edge.neighbors.first);
	const ComputationalCell& left_cell = cells[left_index];
	if(left_cell.stickers.find("dummy")->second)
	  return Conserved();
	const Primitive left = convert_to_primitive(left_cell, eos);
	const Vector2D p = Parallel(edge);
	//	const Primitive right = reflect(left,p);
	const Primitive right = left;
	const Vector2D n = remove_parallel_component
	  (edge.vertices.second - 
	   tess.GetMeshPoint(edge.neighbors.first), p);
	return rotate_solve_rotate_back
	  (rs_, left, right, 0, n, p);
      }
      const size_t left_index =
	static_cast<size_t>(edge.neighbors.first);
      const size_t right_index =
	static_cast<size_t>(edge.neighbors.second);
      const ComputationalCell& left_cell = cells[left_index];
      const ComputationalCell& right_cell = cells[right_index];
      if(left_cell.stickers.find("dummy")->second && 
	 right_cell.stickers.find("dummy")->second)
	return Conserved();
      const Vector2D p = Parallel(edge);
      const Vector2D n = 
	tess.GetMeshPoint(edge.neighbors.second) - 
	tess.GetMeshPoint(edge.neighbors.first);
      const double velocity = Projection
	(tess.CalcFaceVelocity
	 (point_velocities[left_index],
	  point_velocities[right_index],
	  tess.GetCellCM(edge.neighbors.first),
	  tess.GetCellCM(edge.neighbors.second),
	  calc_centroid(edge)),n);			   
      if(left_cell.stickers.find("dummy")->second){
	const Primitive right = 
	  convert_to_primitive(right_cell, eos);
	ComputationalCell ghost;
	ghost.density = right.Density/100;
	ghost.pressure = right.Pressure/100;
	ghost.velocity = Vector2D(0,0);
	const Primitive left = convert_to_primitive(ghost,eos);
	  /*
	  ScalarProd(n,right.Velocity) < 0 ? right : 
	  reflect(right,p);
	  */
	return rotate_solve_rotate_back
	  (rs_,left,right,velocity,n,p);
      }
      if(right_cell.stickers.find("dummy")->second){
	const Primitive left = 
	  convert_to_primitive(left_cell, eos);
	ComputationalCell ghost;
	ghost.density = left.Density/100;
	ghost.pressure = left.Pressure/100;
	ghost.velocity = Vector2D(0,0);
	const Primitive right = convert_to_primitive(ghost,eos);
	  /*
	  ScalarProd(n,left.Velocity)>0 ?
	  left : reflect(left,p);
	  */
	return rotate_solve_rotate_back
	  (rs_,left,right,velocity,n,p);
      }
      const Primitive left = 
	convert_to_primitive(left_cell, eos);
      const Primitive right =
	convert_to_primitive(right_cell, eos);
      return rotate_solve_rotate_back
	(rs_,left,right,velocity,n,p);
    }
  };

  class TimeStepCalculator: public LazyList<double>
  {
  public:

    TimeStepCalculator(const Tessellation& tess,
		       const vector<ComputationalCell>& cells,
		       const EquationOfState& eos,
		       const vector<Vector2D>& point_velocities,
		       const double gm):
      tess_(tess), cells_(cells), 
      point_velocities_(point_velocities), eos_(eos), gm_(gm) {}

    size_t size(void) const 
    {
      return cells_.size();
    }

    double operator[](size_t i) const
    {
      const double kepler_velocity = sqrt(gm_/abs(tess_.GetMeshPoint(static_cast<int>(i))));
      return tess_.GetWidth(static_cast<int>(i))/
	(eos_.dp2c(cells_[i].density, cells_[i].pressure)+
	 abs(cells_[i].velocity)+
	 kepler_velocity+
	 abs(point_velocities_[i]));
    }

  private:
    const Tessellation& tess_;
    const vector<ComputationalCell>& cells_;
    const vector<Vector2D>& point_velocities_;
    const EquationOfState& eos_;
    const double gm_;
  };

  class LogCFL: public TimeStepFunction
  {
  public:

    LogCFL(double cfl, const string& fname, double gm):
      cfl_(cfl), fname_(fname), gm_(gm) {}

    double operator()
    (const Tessellation& tess,
     const vector<ComputationalCell>& cells,
     const EquationOfState& eos,
     const vector<Vector2D>& point_velocities,
     const double /*time*/) const
    {
      const TimeStepCalculator tsc(tess,cells,eos,point_velocities,gm_);
      double res = 10000;
      size_t argmin = 0;
      for(size_t i=0;i<tsc.size();++i){
	const double candidate = tsc[i];
	if(candidate<res && !cells[i].stickers.find("dummy")->second){
	  res = candidate;
	  argmin = i;
	}	
      }
      ofstream f(fname_.c_str());
      f << argmin << endl;
      f << tess.GetMeshPoint(static_cast<int>(argmin)).x << ", ";
      f << tess.GetMeshPoint(static_cast<int>(argmin)).y << endl;
      f << tess.GetWidth(static_cast<int>(argmin)) << endl;
      f << eos.dp2c(cells[argmin].density, cells[argmin].pressure) << endl;
      f << cells[argmin].velocity.x << ", ";
      f << cells[argmin].velocity.y << endl;
      f << point_velocities[argmin].x << ", ";
      f << point_velocities[argmin].y << endl;
      f.close();
      return cfl_*res;
    }

  private:
    const double cfl_;
    const string fname_;
    const double gm_;
  };

  class EntropyFix: public CellUpdater
  {
  public:

    EntropyFix(const double thres=1e-9):
      thres_(thres) {}

    vector<ComputationalCell> operator()
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

  private:
    const double thres_;
  };

  class EdgeLengthCalculator: public LazyList<double>
  {
  public:

    EdgeLengthCalculator(const Tessellation& tess,
			 const PhysicalGeometry& pg):
      tess_(tess), pg_(pg) {}

    size_t size(void) const
    {
      return tess_.getAllEdges().size();
    }

    double operator[](size_t i) const
    {
      return pg_.calcArea(tess_.getAllEdges()[i]);
    }

  private:
    const Tessellation& tess_;
    const PhysicalGeometry& pg_;
  };

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
