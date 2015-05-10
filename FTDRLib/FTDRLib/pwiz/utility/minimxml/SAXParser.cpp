//
// $Id: SAXParser.cpp 2823 2011-06-29 18:37:36Z chambm $
//
//
// Original author: Darren Kessner <darren@proteowizard.org>
//
// Copyright 2007 Spielberg Family Center for Applied Proteomics
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

#define PWIZ_SOURCE

#include "SAXParser.hpp"
#include "pwiz/utility/misc/Std.hpp"
#include "pwiz/utility/misc/random_access_compressed_ifstream.hpp"
#include "boost/regex.hpp"


const string CDATA_begin("![CDATA["), CDATA_end("]]");


//#define PWIZ_USE_BOOST_REGEX 
#ifdef PWIZ_USE_BOOST_REGEX 
#define BOOST_REGEX_MATCH_EXTRA // use boost::regex with multiple matches for tag parsing
#include "boost/regex.hpp"
//
// must link with boost::regex with BOOST_REGEX_MATCH_EXTRA
// defined in boost/regex/user.hpp
//
// boost::regex quick ref:
//   \w (alphanumeric)
//   \s (whitespace)
//   (...) marks a capture
//   (?:...) marks a group with no capture
//   what[0] stores the full matched expression
//   what[n] stores capture n
//
#endif // PWIZ_USE_BOOST_REGEX


namespace pwiz {
namespace minimxml {
namespace SAXParser {
    

namespace {


size_t stripws(string& text)
{
    static const char* ws = "\n\r\t ";

    string::size_type first = text.find_first_not_of(ws);
    text.erase(0, first);

    string::size_type second = text.find_last_not_of(ws);
    if (second+1 < text.size()) text.erase(second+1);

    return first; // return the whitespace stripped from beginning
}


#ifdef PWIZ_USE_BOOST_REGEX

struct ProcessingInstruction 
{
    string name;
    string value; 

    ProcessingInstruction(const string& buffer)
    {
        // ?name value?
        static const boost::regex e("\\?(\\w+)\\s+([^?]+)\\?");

        boost::smatch what; 
        if (!regex_match(buffer,what,e,boost::match_extra) || what.size()!=3)
            throw runtime_error(("[SAXParser::ProcessingInstruction()] Error parsing tag: " + buffer).c_str());

        name = what[1];
        value = what[2];
    }
};

struct StartTag
{
    string name;
    Handler::Attributes attributes;
    bool end;

    StartTag(const string& buffer)
    {
        // tagName (name="value")* '/'?
        static const boost::regex e("([\\w:]+)(?:\\s*([\\w:]+)=\"([^\"]*)\")*\\s*(/)?");

        boost::smatch what; 
        if (!regex_match(buffer,what,e,boost::match_extra) || what.size()!=5)
            throw runtime_error(("[SAXParser::StartTag()] Error parsing tag: " + buffer).c_str());

        if (what.captures(2).size() != what.captures(3).size())
            throw runtime_error("[SAXParser::StartTag()] Error parsing attributes.");

        name = what[1];

        for (unsigned int i=0; i<what.captures(2).size(); i++)
            attributes[what.captures(2)[i]] = what.captures(3)[i];

        end = (what[4]=='/'); 
    }
};

#else // parsing with std library

struct ProcessingInstruction 
{
    string name;
    string value; 

    ProcessingInstruction(const string& buffer)
    {
        istringstream iss(buffer);
        char questionMark = '\0';
        iss >> questionMark >> name;
        if (questionMark != '?') throw runtime_error("[SAXParser::ProcessingInstruction] Error.");
        getline(iss, value, '?');
        stripws(value);
    }
};


string& unescapeXML(string& str)
{
    for (size_t i=0, end=str.size(); i < end; ++i)
    {
        if (str[i] != '&')
            continue;

        // there must be at least three characters after '&' (&lt; or &gt;)
        if (i+3 >= end)
            throw runtime_error("[SAXParser::unescapeXML] Invalid escape sequence.");

        if (str[i+1] == 'l' && str[i+2] == 't' && str[i+3] == ';')
        {
            str[i] = '<';
            str.erase(i+1, 3);
            end -= 3;
        }
        else if (str[i+1] == 'g' && str[i+2] == 't' && str[i+3] == ';')
        {
            str[i] = '>';
            str.erase(i+1, 3);
            end -= 3;
        }
        else if (i+4 < end && str[i+1] == 'a' && str[i+2] == 'm' && str[i+3] == 'p' && str[i+4] == ';')
        {
            str[i] = '&';
            str.erase(i+1, 4);
            end -= 4;
        }
        else if (i+5 < end && str[i+1] == 'q' && str[i+2] == 'u' && str[i+3] == 'o' && str[i+4] == 't' && str[i+5] == ';')
        {
            str[i] = '"';
            str.erase(i+1, 5);
            end -= 5;
        }
        else if (i+5 < end && str[i+1] == 'a' && str[i+2] == 'p' && str[i+3] == 'o' && str[i+4] == 's' && str[i+5] == ';')
        {
            str[i] = '\'';
            str.erase(i+1, 5);
            end -= 5;
        }
        else
            throw runtime_error("[SAXParser::unescapeXML] Invalid escape sequence.");
    }
    return str;
}


const string whitespace_ = " \t\n\r";
const string quote_ = "\"\'";


void parseAttribute(const string& tag, string::size_type& index, Handler::Attributes& attributes, bool unescapeAttributes)
{
    string::size_type indexNameBegin = tag.find_first_not_of(whitespace_, index);
    string::size_type indexNameEnd = tag.find_first_of(whitespace_ + '=', indexNameBegin+1);
    string::size_type indexEquals = tag.find_first_of('=', indexNameEnd);
    string::size_type indexQuoteOpen = tag.find_first_of(quote_, indexEquals+1);
    char quoteChar = tag[indexQuoteOpen];
    string::size_type indexQuoteClose = tag.find_first_of(quoteChar, indexQuoteOpen+1);

    if (indexNameBegin == string::npos ||
        indexNameEnd == string::npos ||
        indexEquals == string::npos ||
        indexQuoteOpen == string::npos ||
        indexQuoteClose == string::npos)
        throw runtime_error("[SAXParser::parseAttribute()] Error at index " 
                            + lexical_cast<string>(index) + ":\n" + tag);

    string name = tag.substr(indexNameBegin, indexNameEnd-indexNameBegin);
    string value = tag.substr(indexQuoteOpen+1, indexQuoteClose-indexQuoteOpen-1);

    if (unescapeAttributes)
        unescapeXML(value);
    attributes[name] = value;
    index = tag.find_first_not_of(whitespace_, indexQuoteClose+1);
}


struct StartTag
{
    string name;
    Handler::Attributes attributes;
    bool end;

    StartTag(const string& buffer, bool unescapeAttributes)
    :   end(false)
    {
        if (buffer[buffer.size()-1] == '/') end = true;

        string::size_type indexNameBegin = buffer.find_first_not_of(whitespace_);
        if (indexNameBegin == string::npos)
            throw runtime_error("[SAXParser::StartTag] Empty buffer.");

        string::size_type indexNameEnd = buffer.find_first_of(whitespace_ + "/", indexNameBegin+1);
        if (indexNameEnd == string::npos)
        {
            name = buffer.substr(indexNameBegin);
            return;
        }

        name = buffer.substr(indexNameBegin, indexNameEnd-indexNameBegin);

        string::size_type index = indexNameEnd + 1;
        string::size_type indexEnd = end ? buffer.size()-1 : buffer.size();
        while (index < indexEnd)
            parseAttribute(buffer, index, attributes, unescapeAttributes);
    }
};

#endif // PWIZ_USE_BOOST_REGEX


struct HandlerInfo
{
    Handler& handler;
    stack<string> names;

    HandlerInfo(Handler& _handler)
    :   handler(_handler)
    {}
};


// HandlerWrangler responsibilities:
// - maintain a Handler stack
// - validate return Status from Handler calls
// - validate element start/end tag matching
class HandlerWrangler : public SAXParser::Handler
{
    public:

    HandlerWrangler(Handler& root)
    {
        handlers_.push(root);
    }

    void verifyNoDelegate(const Status& status)
    {
        if (status.flag==Status::Delegate || status.delegate)
            throw runtime_error("[SAXParser] Illegal return of Status::Delegate.");
    }

    virtual Status processingInstruction(const string& name,
                                         const string& data,
                                         stream_offset position)
    {
        Status status = handlers_.top().handler.processingInstruction(name, data, position);
        verifyNoDelegate(status); 
        return status;
    }

    virtual Status startElement(const string& name,
                                const Attributes& attributes,
                                stream_offset position)
    {
        HandlerInfo& top = handlers_.top();
         
        // element start/end validation

        top.names.push(name);

        // call handler

        Handler::Status status = top.handler.startElement(name, attributes, position);
        if (status.flag != Handler::Status::Delegate)
            return status; 

        // Status::Delegate: let status.delegate handle this message 

        if (!status.delegate) throw runtime_error("[SAXParser] Null delegate.");
        top.names.pop();            
        handlers_.push(*status.delegate);
        return startElement(name, attributes, position);
    }

    virtual Status endElement(const string& name, stream_offset position)
    {
        HandlerInfo& top = handlers_.top();

        // element start/end validation

        if (top.names.empty() || top.names.top()!=name) 
            throw runtime_error(("[SAXParser::ParserWrangler::elementEnd()] Illegal end tag: " + name).c_str()); 

        top.names.pop();

        // call handler

        Status status = top.handler.endElement(name, position);
        verifyNoDelegate(status); 

        // delete handler if we're done with it 

        if (top.names.empty())
        {
            handlers_.pop();
            if (handlers_.empty()) return Status::Done;
        }

        return status;
    }

    virtual Status characters(const string& text, stream_offset position)
    {
        Status status = topHandler().characters(text, position);
        verifyNoDelegate(status);
        return status;
    }

    const Handler& topHandler() const {return handlers_.top().handler;}
    Handler& topHandler() {return handlers_.top().handler;}

    private:
    stack<HandlerInfo> handlers_;
};


} // namespace


bool unbalancedQuote(const string& buffer)
{
    size_t quoteCount = 0;
    string::const_iterator pos = buffer.begin(), end = buffer.end();

    for (; pos != end; ++pos)
        if (*pos == '"' || *pos == '\"')
            ++quoteCount;
    
    return ((quoteCount%2)!=0); // need explicit bool operation to quiet some compilers
}

//
// parse() responsibilities: 
// - stream parsing
// - initiation of events
//   - events are routed through a HandlerWrangler
//   - HandlerWrangler handles any XML/Handler validation
// - return on Handler::Status::Done
//
PWIZ_API_DECL void parse(istream& is, Handler& handler)
{
    using boost::iostreams::position_to_offset;

    HandlerWrangler wrangler(handler);
    Handler::stream_offset position = position_to_offset(is.tellg());

    while (is)
    {
        string buffer;

        // read text up to next tag (may be empty)

        if (!getline(is, buffer, '<')) break;

        // position == beginning of characters

        position += stripws(buffer); 

        // TODO: is it possible to detect when Handler::characters() has been overridden?
        const Handler& topHandler = wrangler.topHandler();
        if (!buffer.empty() && topHandler.parseCharacters)
        {
            if (topHandler.autoUnescapeCharacters)
                unescapeXML(buffer);
            Handler::Status status = wrangler.characters(buffer, position);
            if (status.flag == Handler::Status::Done) return;
        }

        // position == beginning of tag

        position = position_to_offset(is.tellg());
        if (position > 0) position--;

        // read tag

        string temp_buffer;
        bool complete = false, inCDATA = false;
        buffer.clear();

        // If in CDATA, fetch more until the section is ended;
        // else check for unbalanced quotes and fetch more until we have
        // the complete tag.
        do
        {
            if (!getline(is, temp_buffer, '>')) break;
            buffer += temp_buffer;

            inCDATA = inCDATA || bal::starts_with(buffer, CDATA_begin);

            if (inCDATA && !bal::ends_with(buffer, CDATA_end) ||
                !inCDATA && unbalancedQuote(buffer))
                buffer += ">";
            else
                complete = true;
        }
        while(!complete);
        
        stripws(buffer);
        if (buffer.empty())
            throw runtime_error("[SAXParser::parse()] Empty tag."); 

        // switch on tag type

        switch (buffer[0])
        {
            case '?':
            {
                ProcessingInstruction pi(buffer);
                Handler::Status status = wrangler.processingInstruction(pi.name, pi.value, position);
                if (status.flag == Handler::Status::Done) return;
                break; 
            }
            case '/':
            {
                Handler::Status status = wrangler.endElement(buffer.substr(1), position);
                if (status.flag == Handler::Status::Done) return;
                break;
            }
            case '!':
            {
                if (inCDATA)
                {
                    Handler::Status status = wrangler.characters(buffer.substr(CDATA_begin.length(), buffer.size()-CDATA_begin.length()-CDATA_end.length()), position);
                    if (status.flag == Handler::Status::Done) return;
                }
                else if (!bal::starts_with(buffer, "!--") || !bal::ends_with(buffer, "--"))
                    throw runtime_error(("[SAXParser::parse()] Illegal comment: " + buffer).c_str());
                break;
            }
            default: 
            {
                StartTag tag(buffer, handler.autoUnescapeAttributes);

                Handler::Status status = wrangler.startElement(tag.name, tag.attributes, position);
                if (status.flag == Handler::Status::Done) return;
                
                if (tag.end) 
                {
                    status = wrangler.endElement(tag.name, position);
                    if (status.flag == Handler::Status::Done) return;
                }
            }
        }

        // position == after tag end
        position = position_to_offset(is.tellg());
    }
}


} // namespace SAXParser



string xml_root_element(const string& fileheader)
{
    const static boost::regex e("<\\?xml.*?>.*?<([^?!]\\S+?)[\\s>]");

    // convert Unicode to ASCII
    string asciiheader;
    asciiheader.reserve(fileheader.size());
    BOOST_FOREACH(char c, fileheader)
    {
        if(c > 0)
            asciiheader.push_back(c);
    }

    boost::smatch m;
    if (boost::regex_search(asciiheader, m, e))
        return m[1];
    throw runtime_error("[xml_root_element] Root element not found (header is not well-formed XML)");
}

string xml_root_element(istream& is)
{
    char buf[513];
    is.read(buf, 512);
    return xml_root_element(buf);
}

string xml_root_element_from_file(const string& filepath)
{
    pwiz::util::random_access_compressed_ifstream file(filepath.c_str());
    if (!file)
        throw runtime_error("[xml_root_element_from_file] Error opening file");
    return xml_root_element(file);
}


namespace { bool isalnum(char& c) {return std::isalnum(c, std::locale::classic());} }


PWIZ_API_DECL string& decode_xml_id(string& str)
{
    std::istringstream parser;
    for (size_t i=0; i < str.length(); ++i)
    {
        size_t found = str.find("_x00");
        if (found != string::npos &&
            found+6 < str.length() &&
            isalnum(str[found+4]) &&
            isalnum(str[found+5]) &&
            str[found+6] == '_')
        {
            parser.clear(); // reset state
            parser.str(str.substr(found+4, 2));
            int value;
            parser >> std::hex >> value;
            char decoded = (char) value;
            str.replace(found, 7, &decoded, 1);
        }
        else
            break;
    }

    return str;
}


PWIZ_API_DECL string decode_xml_id_copy(const string& str)
{
    string copy(str);
    return decode_xml_id(copy);
}


} // namespace minimxml
} // namespace pwiz


