/*
 * Copyright 2015 Universidad Complutense de Madrid
 *
 * This file is part of Numina
 *
 * Numina is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Numina is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Numina.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

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
