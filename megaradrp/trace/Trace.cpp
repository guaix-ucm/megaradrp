#include <vector>
#include <cstddef>
#include <algorithm>

#include "Trace.h"

#include "fitter.h"

namespace Numina {

  Trace::Trace() 
  {}
  
  void Trace::push_back(double x, double y) {
    xtrace.push_back(x);
    ytrace.push_back(y);
  }

  void Trace::reverse() {
    std::reverse(xtrace.begin(), xtrace.end());
    std::reverse(ytrace.begin(), ytrace.end());
  }

  double Trace::predict(double x) const {

    size_t n = std::min<size_t>(5, xtrace.size());
    Numina::LinearFit mm = Numina::linear_fitter(xtrace.end() - n, xtrace.end(), ytrace.end() - n, ytrace.end());
    return mm.slope * x + mm.intercept;
  }

} // namespace numina
