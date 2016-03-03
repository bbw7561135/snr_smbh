#include <cmath>
#include "constants.hpp"
#include "source/misc/simple_io.hpp"

Constants::Constants(void):
  kilo(1e3),
  centi(1e-2),
  parsec(1),
  meter(3.24e-17*parsec),
  year(1),
  second(year/3.16e7),
  solar_mass(1),
  gram(5.03e-34*solar_mass),
  erg(gram*pow(centi*meter/second,2)),
  newton(kilo*gram*meter/pow(second,2)),
  boltzmann_constant(1.38e-16*erg),
  proton_mass(1.67e-24*gram),
  black_hole_mass(1e7*solar_mass),
  zeta(pow(black_hole_mass/(4.3e6*solar_mass),7./15.)),
  adiabatic_index(5./3.),
  gravitation_constant(6.673e-11*newton*pow(meter/(kilo*gram),2)),
  rho_0(2.2e-22*gram/pow(centi*meter,3)),
  R_0(0.04*parsec*zeta),
  omega_in(1),
  inner_density_prefactor(rho_0*pow(R_0,omega_in)),
  R_b(0.4*parsec*zeta),
  omega_out(3),
  outer_density_prefactor
  (inner_density_prefactor*pow(R_b,omega_out-omega_in)),
  offset(read_number("offset_pc.txt")*parsec),
  wind_speed(700*kilo*meter/second),
  mass_loss_rate
  (read_number
   ("mass_loss_rate_mw.txt")*
   3*3e-3*solar_mass/year/(4*M_PI*pow(0.4*parsec,3))),
  supernova_energy(1e51*erg),
  supernova_radius(0.1*offset),
  supernova_volume((4.*M_PI/3)*pow(supernova_radius,3)),
  supernova_mass(read_number("ejecta_mass_solar_mass.txt")),
  supernova_density(supernova_mass/supernova_volume),
  supernova_pressure((adiabatic_index-1)*supernova_energy/supernova_volume),
  lower_left(parsec*Vector2D(0,-100)),
  upper_right(parsec*Vector2D(100,100)) {}
