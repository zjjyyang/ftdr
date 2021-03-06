//
// $Id: Exception.hpp 2913 2011-08-05 21:20:44Z chambm $
//
//
// Original author: Matt Chambers <matt.chambers .@. vanderbilt.edu>
//
// Copyright 2008 Spielberg Family Center for Applied Proteomics
//   Cedars Sinai Medical Center, Los Angeles, California  90048
// Copyright 2008 Vanderbilt University - Nashville, TN 37232
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

#ifndef _EXCEPTION_HPP_
#define _EXCEPTION_HPP_


#include <stdexcept>


// make debug assertions throw exceptions in MSVC
#ifdef _DEBUG
#include <crtdbg.h>
#include <iostream>
#include <locale>
#include <sstream>
namespace {

inline std::string narrow(const std::wstring& str)
{
    std::ostringstream oss;
    const std::ctype<char>& ctfacet = std::use_facet< std::ctype<char> >(oss.getloc());
    for (size_t i=0; i < str.size(); ++i)
        oss << ctfacet.narrow(str[i], 0);
    return oss.str();
}

inline int CrtReportHook(int reportType, char *message, int *returnValue)
{
    throw std::runtime_error(message);
}

inline int CrtReportHookW(int reportType, wchar_t *message, int *returnValue)
{
    throw std::runtime_error(narrow(message));
}

} // namespace

struct ReportHooker
{
    ReportHooker()
    {
        _CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
        _CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDERR);
        _CrtSetReportHook2(_CRT_RPTHOOK_INSTALL, &CrtReportHook);
        _CrtSetReportHookW2(_CRT_RPTHOOK_INSTALL, &CrtReportHookW);
    }

    ~ReportHooker()
    {
        _CrtSetReportHook2(_CRT_RPTHOOK_REMOVE, &CrtReportHook);
        _CrtSetReportHookW2(_CRT_RPTHOOK_REMOVE, &CrtReportHookW);
    }
};
static ReportHooker reportHooker;
#endif // _DEBUG


// make Boost assertions throw exceptions
#if !defined(NDEBUG)
#define BOOST_ENABLE_ASSERT_HANDLER
#include <sstream>
namespace boost
{
    inline void assertion_failed(char const * expr, char const * function, char const * file, long line) // user defined
    {
        std::ostringstream oss;
        oss << "[" << file << ":" << line << "] Assertion failed: " << expr;
        throw std::runtime_error(oss.str());
    }
} // namespace boost
#endif // !defined(NDEBUG)


#endif // _EXCEPTION_HPP_
