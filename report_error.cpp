#include "report_error.hpp"
#include <iostream>

using namespace std;

void report_error(UniversalError const& eo)
{
  cout << eo.GetErrorMessage() << endl;
  cout.precision(14);
  for(size_t i=0;i<eo.GetFields().size();++i)
    cout << eo.GetFields()[i] << " = "
	 << eo.GetValues()[i] << endl;
}
