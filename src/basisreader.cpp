/*
*
*   PURPOSE: Implements basisreader.hpp, a class for reading basis files.
*
*   DATE           AUTHOR           CHANGES
*   ==================================================================
*   30/08/15       Robert Shaw      Original code.
*   04/09/15       Robert Shaw      Now indexes primitives.
*/
 
#include "basisreader.hpp"
#include "ioutil.hpp"
#include <libecpint/ecp.hpp>
#include "error.hpp"
#include "basis.hpp"
#include <algorithm>
#include <iostream>
#include <cstdlib>
 
// Implement class BasisReader
BasisReader::BasisReader(std::map<int, std::string> ns) : names(ns) {
	char const* tmp = std::getenv( "GAMMA_SHARE" );
	if (tmp == NULL)
		libpath = "./share"; 
	else
		libpath  = std::string(tmp);  
	
	obs_list = read_basis_list("basis");	
	jk_list =  read_basis_list("jkfit");
	mp2_list = read_basis_list("mp2fit"); 
	ecp_list = read_basis_list("ecp"); 
}

std::map<std::string, int> BasisReader::read_basis_list(std::string name) {
	
	std::map<std::string, int> basis_list; 
	
	std::string filename = libpath + "/gamma/basislibrary/" + name + ".list"; 
	pugi::xml_document doc; 
	pugi::xml_parse_result result = doc.load_file(filename.c_str()); 
	pugi::xml_node root = doc.child("root"); 
	
	for (pugi::xml_node item = root.child("Item"); item; item = item.next_sibling("Item")) {
		std::string key = item.attribute("key").value(); 
		std::string value = item.attribute("value").value(); 
		basis_list[key] = std::stoi(value); 
	}
	
	return basis_list; 
	
}

std::vector<libecpint::GaussianShell> BasisReader::read_basis(std::string atom, double* pos, 
												  std::string basis, std::string name, 
												  std::map<std::string, int>& basis_list)
{
	std::vector<libecpint::GaussianShell> basis_set; 
	
	int value; 
	auto it = basis_list.find(basis); 
	if (it != basis_list.end()) {
		value = it->second; 
		std::string filename = 	libpath + "/gamma/basislibrary/bases/build/" + name + std::to_string(value) + ".xml"; 

		pugi::xml_document doc;
		pugi::xml_parse_result result = doc.load_file(filename.c_str());
		pugi::xml_node atom_node = doc.child("root").child(atom.c_str()); 
	
		for (pugi::xml_node shell = atom_node.child("Shell"); shell; shell = shell.next_sibling("Shell")) {

			std::string temp = shell.attribute("lval").value(); 
			int lmult = 0;
			if (temp == "s") { lmult = 0; }
			else if (temp == "sp") { lmult = 1; }
			else if (temp == "p") { lmult = 1; }
			else if (temp == "d") { lmult = 2; }
			else if (temp == "f") { lmult = 3; }
			else if (temp == "g") { lmult = 4; }
			else if (temp == "h") { lmult = 5; }
			else if (temp == "i") { lmult = 6; }
			else if (temp == "k") { lmult = 7; }
			else if (temp == "l") { lmult = 8; }
			
			libecpint::GaussianShell s(pos, lmult); 
			
			for (pugi::xml_node xc = shell.child("xc"); xc; xc = xc.next_sibling("xc")) {
				std::string x = xc.attribute("x").value(); 
				std::string c = xc.attribute("c").value(); 
				s.addPrim(std::stod(x), std::stod(c)); 
			}
			
			basis_set.push_back(s); 
		}
	} 
	
	return basis_set;  
}

void BasisReader::readShellBasis(Basis& b, int q, double *pos, int atom, int type) {
	
	std::string basis_type, atom_name, basis_name;
	std::vector<libecpint::GaussianShell> basis_set;
	
	atom_name = getAtomName(q);  
	
	auto it = names.find(q);
	if (it != names.end()) basis_name = it->second; 
	else {
		it = names.find(0);
		if (it != names.end()) basis_name = it->second;
		else basis_name = "sto-3g";
	}
	std::transform(basis_name.begin(), basis_name.end(), basis_name.begin(), ::tolower);
	
	switch(type) {
		case 1: {
			basis_type = "jkfit"; 
			basis_set = read_basis(atom_name, pos, basis_name, basis_type, jk_list);  
			break;
		} 
		
		case 2: {
			basis_type = "mp2fit"; 
			basis_set = read_basis(atom_name, pos, basis_name, basis_type, mp2_list);
			break; 
		}
		
		default: {
			basis_type = "basis"; 
			basis_set = read_basis(atom_name, pos, basis_name, basis_type, obs_list);
			
			int maxl = 0;
			for (auto& s : basis_set)
				maxl = s.l > maxl ? s.l : maxl; 
			b.setMaxL(maxl); 
		}
	}
	
	if (basis_set.size() > 0) {
		for (auto& shell : basis_set) b.addShell(shell, atom, type); 
	} else {
		throw(Error("BASISLIB", "Could not read basis " + basis_name + " for atom " + atom_name));
	}
}

libecpint::ECP BasisReader::readECP(int q, libecpint::ECPBasis& ecpset, double* center) {

	libecpint::ECP newECP(center); 

	std::string atom_name, basis_name;
	atom_name = getAtomName(q); 
	auto name_it = names.find(-q); 
	if (name_it != names.end()) basis_name = name_it->second; 
	
	auto ecp_it = ecp_list.find(basis_name); 
	if (ecp_it != ecp_list.end()) {
		int value = ecp_it->second; 
		std::string filename = libpath + "/gamma/basislibrary/bases/build/ecp" + std::to_string(value) + ".xml"; 
	
		pugi::xml_document doc;
		pugi::xml_parse_result result = doc.load_file(filename.c_str());
		pugi::xml_node atom_node = doc.child("root").child(atom_name.c_str()); 
		int maxl = std::stoi(atom_node.attribute("maxl").value());
		int ncore = std::stoi(atom_node.attribute("ncore").value()); 
		
		auto it = ecpset.core_electrons.find(q);
		if (it == ecpset.core_electrons.end())
			ecpset.core_electrons[q] = ncore; 
	
		for (pugi::xml_node shell = atom_node.child("Shell"); shell; shell = shell.next_sibling("Shell")) {

			int l = std::stoi(shell.attribute("lval").value());
				
			for (pugi::xml_node nxc = shell.child("nxc"); nxc; nxc = nxc.next_sibling("nxc")) {
				int n = std::stoi(nxc.attribute("n").value()); 
				double x = std::stod(nxc.attribute("x").value()); 
				double c = std::stod(nxc.attribute("c").value()); 
				newECP.addPrimitive(n, l, x, c, true); 
			}
		}
	} else {
		throw(Error("ECPLIB", "Could not read ECP " + basis_name + " for atom " + atom_name)); 
	}
			
	return newECP; 
}
