/*
 *
 *   PURPOSE: Implements ioutil.hpp, a collection of i/o utility functions.
 *
 *   DATE           AUTHOR           CHANGES
 *   ==================================================================
 *   30/08/15       Robert Shaw      Original code.
 *
 */

// Includes
#include "ioutil.hpp"
#include <algorithm>	
#include <vector>
#include <cstdio>
#include <regex>
#include <iostream>
#include <iomanip>
	
std::string names[109] = {"h", "he", "li", "be", "b", "c", "n",
	"o", "f", "ne", "na", "mg", "al", "si", "p", "s", "cl", "ar",
	"k", "ca", "sc", "ti", "v", "cr", "mn", "fe", "co", "ni", "cu",
	"zn", "ga", "ge", "as", "se", "br", "kr", "rb", "sr", "y", "zr",
	"nb", "mo", "tc", "ru", "rh", "pd", "ag", "cd", "in", "sn", "sb",
	"te", "i", "xe", "cs", "ba", "la", "ce", "pr", "nd", "pm", "sm",
	"eu", "gd", "tb", "dy", "ho", "er", "tm", "yb", "lu", "hf", "ta",
	"w", "re", "os", "ir", "pt", "au", "hg", "tl", "pb", "bi", "po",
	"at", "rn", "fr", "ra", "ac", "th", "pa", "u", "np", "pu", "am",
	"cm", "bk", "cf", "es", "fm", "md", "no", "lr", "rf", "db", "sg",
	"bh", "hs", "mt" };
	
const int valences[109] = { 1, 2, 1, 2, 3, 4, 5, 6, 7, 8,
		1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 
		9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 1, 2, 3, 
		4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 
		18, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 
		14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
		26, 27, 28, 29, 30, 31, 32, 1, 2, 3, 4, 5, 6, 7,
		8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
		21, 22, 23}; 
		
const double masses[109] = { 1.0079, 4.0026, 6.941, 9.0122, 10.0811, 
		12.0107, 14.0067, 15.9994, 18.9984, 20.1797, 22.9897,
		24.305, 26.9815, 28.0855, 30.9738, 32.065, 35.453, 39.948,
		39.0983, 40.078, 44.9559, 47.867, 50.9415, 51.9961, 54.938,
	 	55.845, 58.9332, 58.6934, 63.546, 65.39, 69.723, 72.64,
	 	74.9216, 78.96, 79.904, 83.8, 85.4678, 87.62, 88.9059, 91.224,
	 	92.9064, 95.94, 98.0, 101.07, 102.9055, 106.42, 107.8682,
	 	112.411, 114.818, 118.71, 121.76, 127.6, 126.9045, 131.293,
	 	132.9055, 137.327, 138.9055, 140.116, 140.9077, 144.24, 145.0,
		150.36, 151.964, 157.25, 158.9253, 162.5, 164.9303, 167.259,
		168.9342, 173.04, 174.967, 178.49, 180.9479, 183.84, 186.207,
		190.23, 192.217, 195.078, 196.9665, 200.59, 204.3833, 207.2,
		208.9804, 209.0, 210.0, 222.0, 223.0, 226.0, 227.0, 232.0381,
		231.0359, 138.0289, 237.0, 244.0, 243.0, 247.0, 247.0, 251.0,
		252.0, 257.0, 258.0, 259.0, 262.0, 261.0, 262.0, 266.0, 264.0,
		277.0, 268.0 };

// General purpose functions

// Return the atomic mass of an atom with atomic number q
double getAtomMass(int q)
{
	// Array of masses in atomic units for Hydrogen to Meitnerium
	// Ds and all the unun- types don't have reliable masses! 
	
	double mass = q < 1 ? 0 : masses[q-1]; 
 	return mass;
}

int getAtomValence(int q) {
	// Array of ground state valencies for Hydrogen to Meitnerium
		
	int valence = q < 1 ? 0 : valences[q-1];
	return valence;  
}


// Return the text version of atom with atomic number q -
// all caps is used for ease of parsing
// e.g., 6 -> C,   20 -> CA, etc.
std::string getAtomName(int q)
{
	std::string name = q < 1 ? "X" : names[q-1];  
	return name; 
}

// Return the atomic number of an atom with text n
// the reverse of getAtomName
int getAtomCharge(const std::string& n)
{
	// Make it upper case
	std::string name = n;
	std::transform(name.begin(), name.end(), name.begin(), ::tolower);
	
	int q = 1; 
	if (name != "X") { 
	  bool found = false;
	  while (!found && q < 109 ){
		  if(getAtomName(q) == name){
		  	found = true;
		  }
		  q++;
	  }
	}
    return q-1;
}

// Get the text name for a shell of angular momentum l
// e.g. l = 0 -> s,  l = 3 -> f
std::string getShellName(int l)
{
	std::string shell;
	switch(l) {
		case 0: {
			shell = "s";
			break;
		}
		case 1: {
			shell = "p";
			break;
		}	
		case 2: {
			shell = "d";
			break;
		}
		case 3: {
			shell = "f";
			break;
		}
		case 4: {
			shell = "g";
			break;
		}
		case 5: {
			shell = "h";
			break;
		}
		default: {
			shell = "N";
			break;
		}
	}		
	return shell;
}	

bool open_outfile(std::string name, std::ofstream &file) {
	
	bool success = false; 
	
	std::ifstream checkfile(name);
	if (checkfile.is_open()) {
		checkfile.close(); 
		
		int number = 1; 
		std::string new_name; 
		while ( !success && number < 100) {
			new_name = name + "_" + std::to_string(number); 
			checkfile.open(new_name);
			
			if (checkfile.is_open()) {
				checkfile.close();
				number++; 
			} else {
				int chk = std::rename(name.c_str(), new_name.c_str()); 
				success = (chk == 0);
			}
		}
	} 
	
	file.open(name);
	success = file.is_open();
	
	return success; 
}

std::vector<std::string> read_xyz(std::string &filename) {
	std::vector<std::string> geometry; 
	
	std::ifstream xyzfile(filename);
	if (xyzfile.is_open()) {
		std::string line, token;
		size_t pos;
		
		std::regex ws_re("\\s+"); // whitespace
		std::string sep = ", "; 
		
		while(std::getline(xyzfile, line)) {
			    std::vector<std::string> result{
			        std::sregex_token_iterator(line.begin(), line.end(), ws_re, -1), {}
			    };
				if (result.size() == 4) { 
					std::string newline = ""; 
					for (auto& s : result)
						newline += s + sep; 
					geometry.push_back(newline); 
				}
		}
		
		xyzfile.close();
	} 
	
	return geometry;
}

std::vector<int> bounds(int parts, int mem) {
    std::vector<int>bnd;
    int delta = mem / parts;
    int reminder = mem % parts;
    int N1 = 0, N2 = 0;
    bnd.push_back(N1);
    for (int i = 0; i < parts; ++i) {
        N2 = N1 + delta;
        if (i == parts - 1)
            N2 += reminder;
        bnd.push_back(N2);
        N1 = N2;
    }
    return bnd;
}


