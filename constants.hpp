#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP 1

#include "source/tessellation/geometry.hpp"

class Constants
{
public:

  Constants(void);

  // Metric prefix
  const double kilo;
  const double centi;

  // Length / distance
  const double parsec;
  const double meter;

  // Time
  const double year;
  const double second;

  // Mass
  const double solar_mass;
  const double gram;

  // Energy
  const double erg;

  // Force
  const double newton;

  // Parameters
  const double boltzmann_constant;
  const double proton_mass;
  const double black_hole_mass;
  const double zeta;
  const double adiabatic_index;
  const double gravitation_constant;
  const double rho_0;
  const double R_0;
  const double omega_in;
  const double inner_density_prefactor;
  const double R_b;
  const double omega_out;
  const double outer_density_prefactor;
  const double offset;
  const double wind_speed;
  const double supernova_energy;
  const double supernova_radius;
  const double supernova_volume;
  const double supernova_mass;
  const double supernova_density;
  const double supernova_pressure;
  const Vector2D lower_left;
  const Vector2D upper_right;
};

#endif // CONSTANTS_HPP
