/* 
 * 	Copyright (c) 2017 Robert Shaw
 *
 * 	Permission is hereby granted, free of charge, to any person obtaining
 *	a copy of this software and associated documentation files (the
 * 	"Software"), to deal in the Software without restriction, including
 * 	without limitation the rights to use, copy, modify, merge, publish,
 * 	distribute, sublicense, and/or sell copies of the Software, and to
 * 	permit persons to whom the Software is furnished to do so, subject to
 *	the following conditions:
 *
 *	The above copyright notice and this permission notice shall be
 * 	included in all copies or substantial portions of the Software.
 *
 *	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *	EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *	MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *	NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *	LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *	OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *	WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef ERRORHEADERDEF
#define ERRORHEADERDEF

#include <string>
#include <map>

/*! \class Error
	\brief A generic error message.
 */
class Error
{
private:
  std::string code; ///< Shorthand code for error
  std::string msg;  ///< Detailed error message
public:
  Error(); ///< Creates an empty, unspecified error.
  
  /*! Creates a new Error. 
  	  @param c - the short code for the error.	
  	  @param m - the full length error message.
  	*/
  Error(std::string c, std::string m); 
  
  /*! Creates a new Error by looking up the short code against a map
   	  of possible error messages. 
  
  	  @param c - the short code for the error.
  	  @param errormap - a map from codes to full error messages.
   */
  Error(std::string c, std::map<std::string, std::string> errormap);
  
  /*! Prints the error message to the primary output stream (for debugging purposes).
   	  @param full - if True, the full error message is printed; otherwise, just the short code is printed.
   */
  void print(bool full = true); 
  
  std::string getCode() const { return code; } ///< @return The short code for the error.
  std::string getMsg() const { return msg; } ///< @return The full error message.
};

#endif
