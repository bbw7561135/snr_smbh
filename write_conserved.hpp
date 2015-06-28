#ifndef WRITE_CONSERVED_HPP
#define WRITE_CONSERVED_HPP 1

#include "source/newtonian/test_2d/main_loop_2d.hpp"

class WriteConserved: public DiagnosticFunction
{
public:

  WriteConserved(const string& fname);

  void operator()(const hdsim& sim);

  ~WriteConserved(void);

private:
  mutable vector<pair<double,Extensive> > time_cons_;
  const string fname_;
};

#endif // WRITE_CONSERVED_HPP
