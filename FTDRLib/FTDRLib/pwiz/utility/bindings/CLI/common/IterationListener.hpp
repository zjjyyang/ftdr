//
// $Id: IterationListener.hpp 2658 2011-04-22 23:40:14Z chambm $
//
//
// Original author: Matt Chambers <matt.chambers .@. vanderbilt.edu>
//
// Copyright 2011 Vanderbilt University - Nashville, TN 37232
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

#ifndef _ITERATIONLISTENER_HPP_CLI_
#define _ITERATIONLISTENER_HPP_CLI_

#pragma warning( push )
#pragma warning( disable : 4634 4635 )
#include "SharedCLI.hpp"
#include "pwiz/utility/misc/IterationListener.hpp"
#pragma warning( pop )


namespace pwiz {
namespace CLI {
namespace util {


public ref struct IterationListener
{
    enum class Status {Ok, Cancel};

    ref struct UpdateMessage
    {
        UpdateMessage(int iterationIndex, int iterationCount, System::String^ message);

        property System::String^ message;
        property int iterationIndex;
        property int iterationCount;
    };

    IterationListener();

    virtual Status update(UpdateMessage^ updateMessage) {return Status::Ok;}
};

public delegate IterationListener::Status IterationListenerUpdate(IterationListener::UpdateMessage^ updateMessage);


public ref class IterationListenerRegistry
{
    DEFINE_INTERNAL_BASE_CODE(IterationListenerRegistry, pwiz::util::IterationListenerRegistry);
    System::Collections::Generic::Dictionary<IterationListener^,
                                             System::Collections::Generic::KeyValuePair<IterationListenerUpdate^,
                                                                                        System::IntPtr> >^ _listeners;

    public:

    IterationListenerRegistry();

    void addListener(IterationListener^ listener, System::UInt32 iterationPeriod);
    void addListenerWithTimer(IterationListener^ listener, double timePeriod); // seconds
    void removeListener(IterationListener^ listener);

    IterationListener::Status broadcastUpdateMessage(IterationListener::UpdateMessage^ updateMessage);
};


} // namespace util
} // namespace CLI
} // namespace pwiz


#endif // _ITERATIONLISTENER_HPP_CLI_
