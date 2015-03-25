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

#include <vector>
#include <cstddef>
#include <algorithm>

#include "Trace.h"

#include "fitter.h"

namespace Numina {

  Trace::Trace() 
  {}
  
  void Trace::push_back(double x, double y, double p) {
    xtrace.push_back(x);
    ytrace.push_back(y);
    ptrace.push_back(p);
  }

  void Trace::reverse() {
    std::reverse(xtrace.begin(), xtrace.end());
    std::reverse(ytrace.begin(), ytrace.end());
    std::reverse(ptrace.begin(), ptrace.end());
  }

  double Trace::predict(double x) const {

    size_t n = std::min<size_t>(5, xtrace.size());
    Numina::LinearFit mm = Numina::linear_fitter(xtrace.end() - n, xtrace.end(), ytrace.end() - n, ytrace.end());
    return mm.slope * x + mm.intercept;
  }

} // namespace numina
