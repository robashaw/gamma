/* Implementation for error.hpp  
 *
 *     DATE             AUTHOR                CHANGES                             
 *   =======================================================================          
 *     14/08/15         Robert Shaw           Original code             
 *
 */

#include "error.hpp"
#include <iostream>

// Default constructor                                                            
Error::Error()
{
  code = "GEN"; // General error                                                 
  msg = "An unspecified error has occurred.";
}

// Other constructors                                                        
Error::Error(std::string c, std::string m)
{
  code = c;
  msg = m;
}

Error::Error(std::string c, std::map<std::string, std::string> errormap)
{
  code = c;
  // Try and find the code in the errormap                                     
  std::map<std::string, std::string>::iterator it = errormap.find(c);
  if (it != errormap.end()) { // If it exists                                   
    msg = errormap[c]; // Copy in the message                                  
  } else {
    msg = "An unspecified error has occurred."; // Otherwise, generic message   
  }
}

// Print                                                      
void Error::print(bool full)
{
  // Print out code                                                    
  std::cout << code;
  // and message, if required                                                 
  if(full){
    std::cout << ": " << msg << "\n";
  }
}
