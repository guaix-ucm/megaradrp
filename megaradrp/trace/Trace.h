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

#ifndef NU_TRACE_H
#define NU_TRACE_H


#include <vector>

namespace Numina {

class Trace {
  public:
    Trace();
    double predict(double x) const;
    void push_back(double x, double y, double p);
    void reverse();
    std::vector<double> xtrace;
    std::vector<double> ytrace;
    std::vector<double> ptrace;
  };

} // namespace numina

#endif // NU_TRACE_H
