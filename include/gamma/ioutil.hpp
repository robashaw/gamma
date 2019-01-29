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

#ifndef IOUTILHEADERDEF
#define IOUTILHEADERDEF

// General routines
#include <string>
#include <thread>
#include <vector>
#include <fstream>

// Get the atomic mass/name/valency of an atom with atomic number q
double getAtomMass(int q);
std::string getAtomName(int q);
int getAtomValence(int q); 

// Get the atomic number of an atom with name n
int getAtomCharge(const std::string& n);

// Split "mem" into "parts", for multithreading
std::vector<int> bounds(int parts, int mem);

// Get the text name for a shell of angular momentum l
std::string getShellName(int l);
  
bool open_outfile(std::string name, std::ofstream &file);
std::vector<std::string> read_xyz(std::string &filename);  
  
template <typename Lambda>
void parallel_do(Lambda& lambda, int nthreads) {
  std::vector<std::thread> threads;
  for (int thread_id = 0; thread_id != nthreads; ++thread_id) {
    if (thread_id != nthreads - 1)
      threads.push_back(std::thread(lambda, thread_id));
    else
      lambda(thread_id);
  }  // threads_id
  for (int thread_id = 0; thread_id < nthreads - 1; ++thread_id)
    threads[thread_id].join();
}

#endif
