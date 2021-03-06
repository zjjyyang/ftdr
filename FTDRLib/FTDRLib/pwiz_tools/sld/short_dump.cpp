//
// $Id: short_dump.cpp 1195 2009-08-14 22:12:04Z chambm $
//
//
// Original author: Darren Kessner <darren@proteowizard.org>
//
// Copyright 2008 Spielberg Family Center for Applied Proteomics
//   Cedars-Sinai Medical Center, Los Angeles, California  90048
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


#include <iostream>
#include <fstream>


using namespace std;


int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cout << "Usage: short_dump filename\n";
        return 1;
    }

    const char* filename = argv[1];
    ifstream is(filename);
    
    while (is)
    {
        // read 2 bytes at a time
        short x = 0;
        is.read((char*)&x, sizeof(short));
        if (!is) break;

        // print characters if printable, /0 if 0, value otherwise
        if (x>=' ' && x<='~')
            cout << (char)x;
        else if (x==0)
            cout << "\\0 <pos " << is.tellg() << " 0x" << hex << is.tellg() << dec << ">\n";
        else
            cout << "<int " << dec << x << " 0x" << hex << x << dec << ">\n"; 
    }

    return 0;
}

