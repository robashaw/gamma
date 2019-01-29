#include "ProgramController.hpp"
#include "error.hpp"
#include "almoscf.hpp"
#include "cc.hpp"
#include "rpa.hpp"
#include "optimiser.hpp"
#include <libint2.hpp>
#include "eigen_wrapper.hpp"
#include <ctf.hpp>

#include <algorithm>

Option::Option(std::string& line) {
	parse(line); 
}

Option::Option(std::string _name, std::string val) : name(_name), _value(val) {} 

Option::Option(const Option& other) {
	name = other.name; 
	_value = other._value; 
} 

void Option::parse(std::string& line) {
	size_t pos = line.find(','); 
	if (pos != std::string::npos) {
		name = line.substr(0, pos); 
		_value = line.substr(pos+1, line.length());
		if (_value == "true" || _value == "yes" || _value == "on") _value = std::to_string(true);
		else if (_value == "false" || _value == "no" || _value == "off") _value = std::to_string(false);
	}
}

Command::Command(std::vector<std::string>& lines, int _molecule_id, std::string _name) : molecule_id(_molecule_id), name(_name) {
	parse(lines); 
}

Command::Command(const Command& other) {
	options = other.options; 
	molecule_id = other.molecule_id; 
	name = other.name; 
}

void Command::parse(std::vector<std::string>& lines) {
	for (auto& line : lines)
		options.push_back(Option(line)); 
}

Construct::Construct(std::vector<std::string>& lines, std::string& _name) : name(_name) {
	parse(lines); 
}

Construct::Construct(const Construct& other) {
	subconstructs = other.subconstructs;
	content = other.content; 
	name = other.name; 
}

void Construct::parse(std::vector<std::string>& lines) {
	size_t pos; 
	std::string token; 
	for (int i = 0; i < lines.size(); i++) {
		pos = lines[i].find('=');
		if (pos != std::string::npos) { // Subconstruct
			token = lines[i].substr(0, pos); 
			
			std::vector<std::string> sublines; 
			pos = lines[++i].find('}'); 
			bool end_found = (pos != std::string::npos);
			int nsubs = 0; 
			
			while(i < lines.size() && !end_found) {
				sublines.push_back(lines[i]);
				pos = lines[i].find('=');
				if (pos != std::string::npos) {
					nsubs++;
					i++;  
				} else {
					pos = lines[i++].find('}');
					if (pos != std::string::npos && nsubs == 0) {
						end_found = true;
						sublines.pop_back();
						i--;
					} else if (pos != std::string::npos) nsubs--; 
				}
			} 
			
			subconstructs.push_back(Construct(sublines, token)); 
		} else content.push_back(lines[i]); 
	}
}

ProgramController::ProgramController(std::ifstream& input, std::ofstream& output, std::ostream& err, std::string& name) : jobname(name), log(*this, output, err) {
	log.init();
	log.flush();
	
	using namespace std::placeholders;
	command_list["hf"] = std::bind(&ProgramController::call_hf, this, _1, _2); 
	command_list["rhf"] = std::bind(&ProgramController::call_rhf, this, _1, _2);  
	command_list["uhf"] = std::bind(&ProgramController::call_uhf, this, _1, _2); 
	command_list["mp2"] = std::bind(&ProgramController::call_mp2, this, _1, _2); 
	command_list["dfmp2"] = std::bind(&ProgramController::call_dfmp2, this, _1, _2); 
	command_list["rpa"] = std::bind(&ProgramController::call_rpa, this, _1, _2); 
	command_list["ccsd"] = std::bind(&ProgramController::call_ccsd, this, _1, _2); 
	command_list["ralmo"] = std::bind(&ProgramController::call_ralmo, this, _1, _2); 
	command_list["ualmo"] = std::bind(&ProgramController::call_ualmo, this, _1, _2); 
	command_list["optg"] = std::bind(&ProgramController::call_optg, this, _1, _2); 
	command_list["optx"] = std::bind(&ProgramController::call_optx, this, _1, _2);
	command_list["nuctest"] = std::bind(&ProgramController::call_nuctest, this, _1, _2);
	
	parse(input); 
	
	if (!is_option_set("memory")) set_option<double>("memory", 100.0);
	if (!is_option_set("nthreads")) set_option<int>("nthreads", 1);
	if (!is_option_set("printeris")) set_option<bool>("printeris", false);
	if (!is_option_set("printorbs")) set_option<bool>("printorbs", false); 
	if (!is_option_set("direct")) set_option<bool>("direct", false);
	if (!is_option_set("intfile")) set_option<std::string>("intfile", "eris.ints"); 
	if (!is_option_set("thrint")) set_option<double>("thrint", 1e-12);
	if (!is_option_set("bprint")) set_option<bool>("bprint", false); 
	if (!is_option_set("withcore")) set_option<bool>("withcore", false); 
	
	log.init_intfile();
	
}

void ProgramController::cleanLine(std::string& line, bool to_lower) {
	size_t pos = line.find('!');
	if (pos != std::string::npos){
		line.erase(pos, line.length());
	}
	
	line.erase(std::remove(line.begin(), line.end(), ' '), line.end()); 
	line.erase(std::remove(line.begin(), line.end(), '\t'), line.end()); 
	
	if (to_lower)
		std::transform(line.begin(), line.end(), line.begin(), ::tolower); 
}

void ProgramController::parse(std::ifstream& input) {
	std::string line, token;
	size_t pos;
	int curr_molecule = -1; 
	
	while(std::getline(input, line)) {
		// Erase any comments
		cleanLine(line);
		
		// Tokenise
		pos = line.find('=');
		if (pos != std::string::npos) { // Construct
			token = line.substr(0, pos);
			std::vector<std::string> sublines; 
			if (std::getline(input, line)) {
				pos = line.find('}');
				bool end_found = (pos != std::string::npos);
				int nsubs = 0;  
			
				while(!end_found) {
					std::string raw = line;
					cleanLine(line);  
					pos = line.find(',');
					if (pos != std::string::npos) {
						token = line.substr(0, pos); 
						if (token == "geomxyz") {
							cleanLine(raw, false);
							line = raw; 
						}
					} 
						
					sublines.push_back(line);
					pos = line.find('=');
					if (pos != std::string::npos) { 
						nsubs++; 
						if (!std::getline(input, line)) end_found = true;
					} else {
						pos = line.find('}');
						if (pos != std::string::npos && nsubs == 0) {
							end_found = true;
							sublines.pop_back(); 
						} else if (pos != std::string::npos) {
							nsubs--;
							if (!std::getline(input, line)) end_found = true;
						}  else if (!std::getline(input, line)) end_found = true; 
					}
				}  
			}
			constructs.push_back(Construct(sublines, token));
			curr_molecule++; 
		} else {
			pos = line.find('('); 
			if (pos != std::string::npos) { // Command
				token = line.substr(0, pos);   
			 
				auto it = command_list.find(token);
				if (it == command_list.end()) {
					Error e("IO", "Could not find command " + token);
					log.error(e);
				} else { // Read in options
					bool options_end = false; 
					std::vector<std::string> lines; 
					pos = line.find(')');
					if (pos == std::string::npos) {
						if(std::getline(input, line)) {
							while(!options_end) {
								pos = line.find(')'); 
								if (pos != std::string::npos) options_end = true; 
								else {
									cleanLine(line); 
									lines.push_back(line);
									if (!std::getline(input, line)) options_end = true; 
								}
							}
						}
					}
				
					commands.push_back(Command(lines, curr_molecule, token)); 
				}
			} else {
			
				pos = line.find(','); 
				if (pos != std::string::npos) // Global option
					global_options.push_back(Option(line)); 
			
			}	
		}
	}
	
	if (get_option<bool>("debug")) {
		std::cout << "OPTIONS: " << std::endl;
		for (auto& op : global_options) std::cout << op.name << " " << op._value << std::endl;
		std::cout << "\nCOMMANDS: " << std::endl;
		for (auto& c : commands) std::cout << c.name << " " << c.molecule_id << std::endl; 
		std::cout << "\nCONSTRUCTS: " << std::endl; 
		for (auto& c : constructs) {
			std::cout << c.name << std::endl; 
			std::cout << "CONTENT" << std::endl;
			for (auto& line : c.content) std::cout << line << std::endl; 
			std::cout << "SUBCONSTRUCTS" << std::endl; 
			for (auto& sc : c.subconstructs) {
				std::cout << sc.name << std::endl;
				std::cout << "SUBCONTENT" << std::endl;
				for (auto& line : sc.content) std::cout << line << std::endl; 
				std::cout << std::endl;
			} 
			std::cout << std::endl;
		}
	}
}

void ProgramController::run() {
	
	libint2::initialize(); 
	Eigen::setNbThreads(get_option<int>("nthreads")); 
	
	focker = nullptr;
	hf = nullptr;
	mp2obj = nullptr;
	ints = nullptr; 
	
	try {
		// Build molecules
		int curr_molecule = 0;
		for (auto& c : constructs) {
			std::shared_ptr<Molecule> m = std::make_shared<Molecule>(shared_from_this(), c);
			m->buildShellBasis();
			m->buildECPBasis();
			m->calcEnuc();
			log.print(*m, true); 
			m->updateBasisPositions(); 
			log.print(m->getBasis(), m, get_option<bool>("bprint")); 
			log.localTime(); 
			log.flush(); 
			
			ints = std::make_shared<IntegralEngine>(m); 
		
			done_hf = false;
			done_transform = false; 
		
			for (auto& cmd : commands) {
				if(cmd.molecule_id != curr_molecule) continue; 
				else {
					auto it = command_list.find(cmd.name);
					if (it != command_list.end())  {
						std::function<void(Command&, SharedMolecule)> f = it->second;
						f(cmd, m); 
					}
				}
			}
		
			curr_molecule++; 
		}
	} catch (Error e) {
		log.error(e); 
	}
	libint2::finalize(); 
	log.finalise(); 
	
}

void ProgramController::call_hf(Command& c, SharedMolecule m) {
	if (m->getMultiplicity() > 1 || m->getNel()%2 != 0)
		call_uhf(c, m);
	else 
		call_rhf(c, m); 
}

void ProgramController::call_rhf(Command& c, SharedMolecule m) {
	if(!c.is_option_set("diis")) c.set_option<bool>("diis", true); 
	if(!c.is_option_set("maxdiis")) c.set_option<int>("maxdiis", 8);
	if(!c.is_option_set("maxiter")) c.set_option<int>("maxiter", 40);
	if(!c.is_option_set("enconverge")) c.set_option<double>("enconverge", 1e-7); 
	if(!c.is_option_set("densconverge")) c.set_option<double>("densconverge", 1e-4);  
	if(!c.is_option_set("precision")) c.set_option<double>("precision", 1e-12); 
	if(!c.is_option_set("momap")) c.set_option<bool>("momap", false);
	if(!c.is_option_set("mapfile")) c.set_option<std::string>("mapfile", "mo.map");
	if(!c.is_option_set("fineness")) c.set_option<int>("fineness", 50); 
	if(!c.is_option_set("orbital")) c.set_option<int>("orbital", 0);
	if(!c.is_option_set("df")) c.set_option<bool>("df", false); 
	if(!c.is_option_set("guess")) c.set_option<std::string>("guess", "soad"); 
	
	focker = std::make_shared<Fock>(c, *ints, m);
	hf = std::make_shared<SCF>(c, m, *focker); 
	hf->rhf(); 
	done_hf = true;
}

void ProgramController::call_uhf(Command& c, SharedMolecule m) {
	if(!c.is_option_set("diis")) c.set_option<bool>("diis", true); 
	if(!c.is_option_set("maxdiis")) c.set_option<int>("maxdiis", 8);
	if(!c.is_option_set("maxiter")) c.set_option<int>("maxiter", 40);
	if(!c.is_option_set("enconverge")) c.set_option<double>("enconverge", 1e-7); 
	if(!c.is_option_set("densconverge")) c.set_option<double>("densconverge", 1e-4);  
	if(!c.is_option_set("precision")) c.set_option<double>("precision", 1e-12); 
	if(!c.is_option_set("df")) c.set_option<bool>("df", false); 
	if(!c.is_option_set("guess")) c.set_option<std::string>("guess", "soad"); 

	focker = std::make_shared<UnrestrictedFock>(c, *ints, m);
	hf = std::make_shared<SCF>(c, m, *focker); 
	hf->uhf(); 
	done_hf = true;
}

void ProgramController::call_mp2(Command& c, SharedMolecule m) {
	if(done_hf) { 
		mp2obj = std::make_shared<MP2>(*focker); 
		
		log.title("MP2 CALCULATION"); 
		mp2obj->tensormp2(); 
		log.print("Integral transformation complete.\n");
		log.localTime();

		log.result("MP2 Energy Correction", mp2obj->getEnergy(), "Hartree");
		log.result("Total Energy = ", hf->getEnergy() + mp2obj->getEnergy(), "Hartree");
		
	} else {
		Error e("MP2", "HF is required before MP2 can be done.");
		log.error(e);
	}
}

void ProgramController::call_dfmp2(Command& c, SharedMolecule m) {
	if(done_hf) { 
		mp2obj = std::make_shared<MP2>(*focker); 
		
		log.title("DF-MP2 CALCULATION"); 
		mp2obj->dfmp2(true); 

		log.result("DF-MP2 Energy Correction", mp2obj->getEnergy(), "Hartree");
		log.result("Total Energy = ", hf->getEnergy() + mp2obj->getEnergy(), "Hartree");
		
	} else {
		Error e("MP2", "HF is required before DF-MP2 can be done.");
		log.error(e);
	}
}

void ProgramController::call_rpa(Command& c, SharedMolecule m) {
	if(!c.is_option_set("sosex")) c.set_option<bool>("sosex", true); 
	if(!c.is_option_set("longrange")) c.set_option<bool>("longrange", false); 
	if(!c.is_option_set("mu")) c.set_option<double>("mu", 0.5); 
	if(!c.is_option_set("iterative")) c.set_option<bool>("iterative", false);
	
	if(done_hf) { 
		CTF::World dw(MPI_COMM_WORLD); 
		RPA rpa(c, *focker, focker->getHCore().rows(), m->getNel()/2, dw); 
		
		log.title("RPA CALCULATION"); 
		rpa.compute();  
		log.localTime();

		log.result("RPA Energy Correction", rpa.getEnergy(), "Hartree");
		log.result("Total Energy = ", hf->getEnergy() + rpa.getEnergy(), "Hartree");
		
	} else {
		Error e("RPA", "HF is required before RPA can be done.");
		log.error(e);
	}
}

void ProgramController::call_ccsd(Command& c, SharedMolecule m) {
	if(!c.is_option_set("diis")) c.set_option<bool>("diis", true);
	if(!c.is_option_set("triples")) c.set_option<bool>("triples", false);
	if(!c.is_option_set("maxdiis")) c.set_option<int>("maxdiis", 5);
	if(!c.is_option_set("maxiter")) c.set_option<int>("maxiter", 30); 
	if(!c.is_option_set("enconverge")) c.set_option<double>("enconverge", 1e-7); 
	if(!c.is_option_set("densconverge")) c.set_option<double>("densconverge", 1e-4); 
	
	if (done_hf) {
		if(!done_transform) {
			mp2obj = std::make_shared<MP2>(*focker); 
			mp2obj->cctrans(); 
			done_transform = true;
		}
		
		CCSD ccobj(c, *mp2obj); 
		ccobj.compute();
		log.result("Total Energy = ", hf->getEnergy() + ccobj.getEnergy() + ccobj.getETriples(), "Hartree");
	} 
}

void ProgramController::call_ralmo(Command& c, SharedMolecule m) {
	if(!c.is_option_set("diis")) c.set_option<bool>("diis", true);
	if(!c.is_option_set("maxdiis")) c.set_option<int>("maxdiis", 6);
	if(!c.is_option_set("enconverge")) c.set_option<double>("enconverge", 1e-6); 
	if(!c.is_option_set("densconverge")) c.set_option<double>("densconverge", 1e-4);  
	if(!c.is_option_set("maxiter")) c.set_option<int>("maxiter", 40);
	if(!c.is_option_set("perturb")) c.set_option<int>("perturb", 1);
	if(!c.is_option_set("precision")) c.set_option<double>("precision", 1e-12);
	if(!c.is_option_set("rpa")) c.set_option<bool>("rpa", false);  
	if(!c.is_option_set("rpax")) c.set_option<int>("rpax", 2); 
	if(!c.is_option_set("longrange")) c.set_option<bool>("longrange", false); 
	if(!c.is_option_set("mu")) c.set_option<double>("mu", 0.5); 
	if(!c.is_option_set("iterative")) c.set_option<bool>("iterative", true);
	if(!c.is_option_set("pairwise")) c.set_option<bool>("pairwise", true); 
	if(!c.is_option_set("local")) c.set_option<bool>("local", true);
	if(!c.is_option_set("xcorrect")) c.set_option<bool>("xcorrect", false); 
	if(!c.is_option_set("rcutoff")) c.set_option<double>("rcutoff", 15.0);
	if(!c.is_option_set("mothresh")) c.set_option<double>("mothresh", 1e-6);
	if(!c.is_option_set("fitthresh")) c.set_option<double>("fitthresh", 0.05);
	if(!c.is_option_set("rthresh")) c.set_option<double>("rthresh", 7.0);
	if(!c.is_option_set("df")) c.set_option<bool>("df", false); 
	if(!c.is_option_set("guess")) c.set_option<std::string>("guess", "soad"); 
	if(!c.is_option_set("dprint")) c.set_option<bool>("dprint", false);

	focker = std::make_shared<Fock>(c, *ints, m); 
	ALMOSCF almo(c, m, *focker); 
	almo.rscf(); 
}

void ProgramController::call_ualmo(Command& c, SharedMolecule m) {
	if(!c.is_option_set("diis")) c.set_option<bool>("diis", true);
	if(!c.is_option_set("maxdiis")) c.set_option<int>("maxdiis", 6);
	if(!c.is_option_set("enconverge")) c.set_option<double>("enconverge", 1e-6); 
	if(!c.is_option_set("densconverge")) c.set_option<double>("densconverge", 1e-4);  
	if(!c.is_option_set("maxiter")) c.set_option<int>("maxiter", 40);
	if(!c.is_option_set("perturb")) c.set_option<int>("perturb", 1);
	if(!c.is_option_set("precision")) c.set_option<double>("precision", 1e-12);
	if(!c.is_option_set("rpa")) c.set_option<bool>("rpa", false);  
	if(!c.is_option_set("rpax")) c.set_option<int>("rpax", 2); 
	if(!c.is_option_set("longrange")) c.set_option<bool>("longrange", false); 
	if(!c.is_option_set("mu")) c.set_option<double>("mu", 0.5); 
	if(!c.is_option_set("iterative")) c.set_option<bool>("iterative", true);
	if(!c.is_option_set("pairwise")) c.set_option<bool>("pairwise", true); 
	if(!c.is_option_set("local")) c.set_option<bool>("local", true);
	if(!c.is_option_set("xcorrect")) c.set_option<bool>("xcorrect", false); 
	if(!c.is_option_set("rcutoff")) c.set_option<double>("rcutoff", 15.0);
	if(!c.is_option_set("mothresh")) c.set_option<double>("mothresh", 1e-6);
	if(!c.is_option_set("fitthresh")) c.set_option<double>("fitthresh", 0.05);
	if(!c.is_option_set("rthresh")) c.set_option<double>("rthresh", 7.0);
	if(!c.is_option_set("df")) c.set_option<bool>("df", false); 
	if(!c.is_option_set("guess")) c.set_option<std::string>("guess", "soad"); 
	if(!c.is_option_set("dprint")) c.set_option<bool>("dprint", false);

	focker = std::make_shared<UnrestrictedFock>(c, *ints, m);
	ALMOSCF almo(c, m, *focker); 
	almo.uscf(); 
}

void ProgramController::call_optg(Command& c, SharedMolecule m) {
	if(!c.is_option_set("diis")) c.set_option<bool>("diis", true); 
	if(!c.is_option_set("maxdiis")) c.set_option<int>("maxdiis", 8);
	if(!c.is_option_set("maxiter")) c.set_option<int>("maxiter", 40);
	if(!c.is_option_set("enconverge")) c.set_option<double>("enconverge", 1e-7); 
	if(!c.is_option_set("densconverge")) c.set_option<double>("densconverge", 1e-4); 
	if(!c.is_option_set("gradconverge")) c.set_option<double>("gradconverge", 1e-4);  
	if(!c.is_option_set("precision")) c.set_option<double>("precision", 1e-12);
	if(!c.is_option_set("trust")) c.set_option<double>("trust", 0.1); 
	if(!c.is_option_set("freq")) c.set_option<bool>("freq", false); 
	if(!c.is_option_set("modes")) c.set_option<bool>("modes", false); 
	if(!c.is_option_set("guess")) c.set_option<std::string>("guess", "soad"); 
	
	RHFOptimiser optim(c, m);
	optim.optimise();  
}

void ProgramController::call_optx(Command& c, SharedMolecule m) {
	if(!c.is_option_set("diis")) c.set_option<bool>("diis", true); 
	if(!c.is_option_set("maxdiis")) c.set_option<int>("maxdiis", 8);
	if(!c.is_option_set("maxiter")) c.set_option<int>("maxiter", 40);
	if(!c.is_option_set("enconverge")) c.set_option<double>("enconverge", 1e-7); 
	if(!c.is_option_set("densconverge")) c.set_option<double>("densconverge", 1e-4); 
	if(!c.is_option_set("gradconverge")) c.set_option<double>("gradconverge", 1e-4);  
	if(!c.is_option_set("precision")) c.set_option<double>("precision", 1e-12);
	if(!c.is_option_set("trust")) c.set_option<double>("trust", 5.0);
	if(!c.is_option_set("active")) c.set_option<std::string>("active", "all"); 
	if(!c.is_option_set("momap")) c.set_option<bool>("momap", false);
	if(!c.is_option_set("mapfile")) c.set_option<std::string>("mapfile", "mo.map");
	if(!c.is_option_set("fineness")) c.set_option<int>("fineness", 50); 
	if(!c.is_option_set("orbital")) c.set_option<int>("orbital", 0);
	if(!c.is_option_set("mp2")) c.set_option<bool>("mp2", false); 
	if(!c.is_option_set("guess")) c.set_option<std::string>("guess", "hcore"); 
	 
	RHFOptimiser optim(c, m); 
	optim.exponents(); 
}

void ProgramController::call_nuctest(Command& c, SharedMolecule m) {
	
	if(!c.is_option_set("enconverge")) c.set_option<double>("enconverge", 1e-7); 
	if(!c.is_option_set("densconverge")) c.set_option<double>("densconverge", 1e-4); 
	if(!c.is_option_set("maxiter")) c.set_option<int>("maxiter", 40);
	
	Atom& a = m->getAtom(0); 
	int Z = a.getCharge(); 
	double mZ = a.getMass() * 1836.15267389;
	 
	
	Matrix TE = ints->getKinetic();
	Matrix TN = ints->getKinetic() / mZ; 
	Matrix S = ints->getOverlap(); 
	
	int nbfs = TE.rows();
	int nocc = m->getNel() / 2; 
	
	EigenSolver es_orthog(S);
	Matrix U = es_orthog.eigenvectors();
	Vector lambda = es_orthog.eigenvalues();  
	Matrix orthog = Matrix::Zero(nbfs, nbfs);
	for (int i = 0; i < nbfs; i++)
		orthog(i, i) = 1.0/(std::sqrt(lambda(i)));
	orthog = U * orthog * U.transpose();
	
	Matrix VE = Matrix::Zero(nbfs, nbfs);
	Matrix VE_tilde = Matrix::Zero(nbfs, nbfs);
	Matrix VN_tilde = Matrix::Zero(nbfs, nbfs); 
	Matrix FE_mo = orthog.transpose() * TE * orthog;
	Matrix FE_ao = Matrix::Zero(nbfs, nbfs); 
	Matrix FN_mo = orthog.transpose() * TN * orthog; 
	Matrix FN_ao = Matrix::Zero(nbfs, nbfs);	
	
	EigenSolver es_E(FE_mo);
	Matrix CP_E = orthog * es_E.eigenvectors();
	Vector eps_E = es_E.eigenvalues();
	EigenSolver es_N(FN_mo);
	Matrix CP_N = orthog * es_N.eigenvectors();
	Vector eps_N = es_N.eigenvalues(); 
	
	Matrix DE = 2.0 * CP_E.block(0, 0, nbfs, nocc) * CP_E.transpose().block(0, 0, nocc, nbfs); 
	Matrix DN = CP_N.block(0, 0, nbfs, 1) * CP_N.transpose().block(0, 0, 1, nbfs); 
	
	for (int u = 0; u < nbfs; u++) {
		for (int v = 0; v < nbfs ; v++) {
			for (int s = 0; s < nbfs; s++) {
				for (int l = 0; l < nbfs; l++) {
					VE(u, v) += 0.5 * DE(l, s)* (ints->getERI(u, v, s, l) - 0.5*ints->getERI(u, l, s, v));
					VN_tilde(u, v) -= Z * DE(l, s) * ints->getERI(u, v, s, l); 
					VE_tilde(u, v) -= Z * DN(l, s) * ints->getERI(s, l, u, v);
				}
			}
		}
	}
	
	FE_ao = TE + 2.0*VE + VE_tilde; 
	FN_ao = TN + VN_tilde; 
	FE_mo = orthog.transpose() * FE_ao * orthog;
	FN_mo = orthog.transpose() * FN_ao * orthog; 
	
	double en_elec = ((TE + VE + VE_tilde) * DE).trace(); 
	double en_kin_nuc = (TN * DN).trace(); 
	double energy = en_elec + en_kin_nuc; 
	
	int MAXITER = c.get_option<int>("maxiter");
	double CONVERGE = c.get_option<double>("enconverge"); 
	
	int iter = 0;
	bool converged = false; 
	
	log.title("Nuclear SCF test");
	log.initIteration(); 
	double old_energy = energy; 
	
	while (!converged && iter < MAXITER) {
		
		es_E.compute(FE_mo);
		CP_E = orthog * es_E.eigenvectors();
		eps_E = es_E.eigenvalues();
		es_N.compute(FN_mo);
		CP_N = orthog * es_N.eigenvectors();
		eps_N = es_N.eigenvalues(); 
	
		DE = 2.0 * CP_E.block(0, 0, nbfs, nocc) * CP_E.transpose().block(0, 0, nocc, nbfs); 
		DN = CP_N.block(0, 0, nbfs, 1) * CP_N.transpose().block(0, 0, 1, nbfs); 
		
		VE = Matrix::Zero(nbfs, nbfs);
		VN_tilde = Matrix::Zero(nbfs, nbfs);
		VE_tilde = Matrix::Zero(nbfs, nbfs);
		for (int u = 0; u < nbfs; u++) {
			for (int v = 0; v < nbfs ; v++) {
				for (int s = 0; s < nbfs; s++) {
					for (int l = 0; l < nbfs; l++) {
						VE(u, v) += 0.5 * DE(l, s)* (ints->getERI(u, v, s, l) - 0.5*ints->getERI(u, l, s, v));
						VN_tilde(u, v) -= Z * DE(l, s) * ints->getERI(u, v, s, l); 
						VE_tilde(u, v) -= Z * DN(l, s) * ints->getERI(s, l, u, v); 
					}
				}
			}
		}
	
		FE_ao = TE + 2.0*VE + VE_tilde; 
		FN_ao = TN + VN_tilde; 
		FE_mo = orthog.transpose() * FE_ao * orthog;
		FN_mo = orthog.transpose() * FN_ao * orthog; 
	
		en_elec = ((TE + VE + VE_tilde) * DE).trace(); 
		en_kin_nuc = (TN * DN).trace(); 
		energy = en_elec + en_kin_nuc; 
		double delta_e = energy - old_energy; 
		old_energy = energy; 
		
		log.iteration(iter, energy, delta_e, en_elec);
		converged = fabs(delta_e) < CONVERGE;
		iter++; 
	}
	
	es_E.compute(FE_mo);
	CP_E = orthog * es_E.eigenvectors();
	eps_E = es_E.eigenvalues();
	es_N.compute(FN_mo);
	CP_N = orthog * es_N.eigenvectors();
	eps_N = es_N.eigenvalues(); 
	
	if (!converged)
		log.result("Failed to converge.");
	else {
		log.print("\nElectronic energy (Hartree) = " + std::to_string(en_elec));
		log.print("\nNuclear kinetic energy (Hartree) = " + std::to_string(en_kin_nuc)); 
		log.print("\n");
		log.orbitals(eps_E, nocc*2, false); 
		log.print(eps_N);
		log.result("Energy", energy, "Hartree"); 
	}
	
}


void ProgramController::runmp2(MP2& mp2obj, SCF& hf, bool calc) {
	if (calc)
		log.title("MP2 CALCULATION");
	else
		log.title("INTEGRAL TRANSFORMATION");
	
	mp2obj.transformIntegrals();
	log.print("Integral transformation complete.\n");
	log.localTime();
	if (calc) {
		mp2obj.calculateEnergy();
		log.result("MP2 Energy Correction", mp2obj.getEnergy(), "Hartree");
		log.result("Total Energy = ", hf.getEnergy() + mp2obj.getEnergy(), "Hartree");
	}
}

ProgramController& ProgramController::operator=(const ProgramController& other) {
	global_options = other.global_options;
	commands = other.commands;
	constructs = other.constructs;
	command_list = other.command_list;
	focker = other.focker;
	hf = other.hf;
	mp2obj = other.mp2obj;
	ints = other.ints;
	done_hf = other.done_hf;
	done_transform = other.done_transform; 
	jobname = other.jobname; 
	log = other.log;
	return *this;
}
