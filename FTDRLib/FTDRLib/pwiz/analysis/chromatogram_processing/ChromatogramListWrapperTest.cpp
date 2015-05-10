//
// $Id: ChromatogramListWrapperTest.cpp 2051 2010-06-15 18:39:13Z chambm $
//
//
// Original author: Eric Purser <Eric.Purser .@. Vanderbilt.edu>
//
// Copyright 2008 Spielberg Family Center for Applied Proteomics
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


#include "ChromatogramListWrapper.hpp"
#include "pwiz/utility/misc/unit.hpp"
#include "pwiz/utility/misc/Std.hpp"

using namespace pwiz::analysis;
using namespace pwiz::cv;
using namespace pwiz::msdata;
using namespace pwiz::util;


class MyWrapper : public ChromatogramListWrapper
{
    public:

    MyWrapper(const ChromatogramListPtr& inner)
    :   ChromatogramListWrapper(inner)
    {}

    void verifySize(size_t size)
    {
        // verify that we can see inner_ 
        unit_assert(size == inner_->size());
    }
};


void test()
{
    ChromatogramListSimplePtr simple(new ChromatogramListSimple);

    const size_t chromatogramCount = 10;
    for (size_t i=0; i<chromatogramCount; i++)
    {
        simple->chromatograms.push_back(ChromatogramPtr(new Chromatogram));
        Chromatogram& s = *simple->chromatograms.back();
        s.index = i;
        s.id = "S" + lexical_cast<string>(i);
        s.nativeID = lexical_cast<string>(i);
    }

    shared_ptr<MyWrapper> wrapper(new MyWrapper(simple)); 

    // make sure we're getting what we expect

    wrapper->verifySize(10);
    unit_assert(wrapper->size() == 10);
    for (size_t i=0; i<chromatogramCount; i++)
    {
        string id = "S" + lexical_cast<string>(i);
        string nativeID = lexical_cast<string>(i);

        unit_assert(wrapper->find(id) == i);
        unit_assert(wrapper->findNative(nativeID) == i);

        const ChromatogramIdentity& identity = wrapper->chromatogramIdentity(i);
        unit_assert(identity.id == id);
        unit_assert(identity.nativeID == nativeID);

        ChromatogramPtr s = wrapper->chromatogram(i);
        unit_assert(s->id == id);
        unit_assert(s->nativeID == nativeID);
    }
}


int main()
{
    try
    {
        test();
        return 0;
    }
    catch (exception& e)
    {
        cerr << e.what() << endl;
        return 1;
    }
    catch (...)
    {
        cerr << "Caught unknown exception.\n";
        return 1;
    }
}


