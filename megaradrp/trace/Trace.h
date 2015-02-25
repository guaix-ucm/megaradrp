#ifndef NU_TRACE_H
#define NU_TRACE_H


#include <vector>

namespace Numina {

class Trace {
  public:
    Trace();
    double predict(double x) const;
    void push_back(double x, double y);
    void reverse();
    std::vector<double> xtrace;
    std::vector<double> ytrace;
  };

} // namespace numina

#endif // NU_TRACE_H
