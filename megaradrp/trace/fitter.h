#ifndef NU_FITTER_H
#define NU_FITTER_H

#include <numeric>
#include <vector>
#include <stdexcept>

namespace Numina {

struct LinearFit {
  double slope;
  double intercept;
  bool fit;
};

template <class Iterator>
LinearFit linear_fitter(Iterator x1, Iterator x2, Iterator y1, Iterator y2) {
   typedef typename std::iterator_traits<Iterator>::difference_type diff_t;
   diff_t xn = std::distance(x1, x2);
   diff_t yn = std::distance(y1, y2);
   
   if (xn != yn)
     throw std::invalid_argument("XN must be == YN");
   
   if (xn < 1)
     throw std::invalid_argument("XN must be > 2");
   
   if (xn == 1) {
     LinearFit result;
     result.slope = 0.0;
     result.intercept = *y1;
     result.fit = false;
     return result;
   }
   double xm = std::accumulate(x1, x2, 0.0) / xn;
   double ym = std::accumulate(y1, y2, 0.0) / yn;

   double acc1 = 0.0;
   double acc2 = 0.0;

   for(diff_t j = 0; j < xn; j++) {
     acc1 += (*(x1 + j) - xm) * (*(y1 + j) - ym);
   }

   for(diff_t j = 0; j < xn; j++) {
     acc2 += (*(x1 + j) - xm) * (*(x1 + j) - ym);
   }

   LinearFit result;
   result.fit = true;
   result.slope = acc1 / acc2;
   result.intercept = ym - result.slope * xm;
   return result;
}


} // namespace Numina

#endif // NU_FITTER_H
