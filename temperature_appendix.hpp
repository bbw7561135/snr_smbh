#ifndef TEMPERATURE_APPENDIX_HPP
#define TEMPERATURE_APPENDIX_HPP 1

#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"

class TemperatureAppendix: public DiagnosticAppendix
{
public:

  TemperatureAppendix(double m2k);

  string getName(void) const;

  vector<double> operator()(const hdsim& sim) const;

private:

  //! \brief Ratio between atomic mass and the Boltzmann constant
  double m2k_;
};

#endif // TEMPERATURE_APPENDIX_HPP
