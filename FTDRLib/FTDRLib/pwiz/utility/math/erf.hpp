//
// $Id: erf.hpp 1195 2009-08-14 22:12:04Z chambm $
//
//
// Original author: Darren Kessner <darren@proteowizard.org>
//
// Copyright 2006 Louis Warschaw Prostate Cancer Center
//   Cedars Sinai Medical Center, Los Angeles, California  90048
//
// Licensed under the Apache License, Version 2.0 (the "License"); 
// you may not use this file except in compliance with the License. 
// You may obtain a copy of the License at 
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software 
// distributed under the License is distributed on an "AS IS" BASIS, 
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
// See the License for the specific language governing permissions and 
// limitations under the License.
//


#ifndef _ERF_HPP_
#define _ERF_HPP_


#include "pwiz/utility/misc/Export.hpp"
#include <complex>


namespace pwiz {
namespace math {


/// real error function; calls gcc-provided erf, complex version (below) on msvc
PWIZ_API_DECL double erf(double x);

/// complex error function 
PWIZ_API_DECL std::complex<double> erf(const std::complex<double>& z);


/// series implementation for testing
std::complex<double> erf_series2(const std::complex<double>& z);


} // namespace math
} // namespace pwiz


#endif // _ERF_HPP_

