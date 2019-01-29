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

#ifndef PROGRAM_CONTROLLER_HEAD
#define PROGRAM_CONTROLLER_HEAD

#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <functional>
#include <sstream>
#include "logger.hpp"
#include "molecule.hpp"
#include "fock.hpp"
#include "scf.hpp"
#include "mp2.hpp"
#include "integrals.hpp"

/*! \struct Option
	\brief	An option for the program or a command.

	Input of the form 
	     name, value
	will be interpreted as an Option. This struct parses and 
	stores these options, accepting any value type that can be
	piped from a stringstream.
 */
struct Option {	
	/// Creates an empty option
	Option() { }
	
	/*! Parses a line of input into an Option.
		@param line - a reference to a string containing the input to be parsed. 
	 */
	Option(std::string& line);
	
	/// Initialises an Option with the given name and value.
	Option(std::string _name, std::string val); 
	
	/// Copy constructor
	Option(const Option& other); 
	
	std::string name; ///< The identifier for this option.
	std::string _value; ///< The type-ambiguous value of the option.
	
	/*! Returns the value of the Option with the specified type, if it can be converted.
		@tparam T - the type of the value; must be pipeable from a stringstream.
		@return The value of the Option as type T.
	 */
	template <typename T>
	T get_value() {
		std::stringstream ss(_value);
		T convertedValue;
		if ( ss >> convertedValue ) return convertedValue;
		else throw std::runtime_error("conversion failed");
	}
	
	/*! Sets the value of the Option with the specified type, if it can be converted.
		@tparam T - the type of the value; must be pipeable into an ostringstream.
	 */
	template <typename T>
	void set_value(T val) {
		std::ostringstream os;
		os << val; 
		_value = os.str(); 
	}
	
	/*! Parses a string of the form "name, value" into an Option.
	 @param line - a reference to the string to be parsed.
	 */
	void parse(std::string& line); 
};

/*! \class Command
	\brief A command for the program, optionally with Options.

	Any input of the form
		name({OPTIONS})
	will be interpreted as a command, and the options parsed into it. 
	
	These commands are then passed to the relevant part of the program. 
*/
class Command {
private:
	std::vector<Option> options; ///< A list of all the options for the command 
	
public:
	/// Creates an empty command named "default" with no molecule attached.
	Command() : molecule_id(-1), name("default") { } 
	
	/*! Create a command with the given name and molecule, and parse the options into it.
		@param lines - strings that can be parsed into Options.
		@param molecule_id - the identifier of the molecule this command is attached to.
		@param name - the name of the command.
	 */
	Command(std::vector<std::string>& lines, int molecule_id, std::string name);
	
	/// Copy constructor
	Command(const Command& other);  
	
	/*! Takes a collection of strings and creates an Option from each.
		@param lines - strings that can be converted into Options.
	 */
	virtual void parse(std::vector<std::string>& lines);
	
	/*! Retrieves the Option with the given name as the type specified.   
		If the Option is not found in the list of options for this command,
		it returns the default value for the type specified (e.g. 0 for ints, etc.).
		@tparam T - the type of return value for the option.
		@param name - the name of the Option wanted.
		@return The value of the Option if found, or the default value for type T if not.
	 */
	template <typename T>
	T get_option(std::string name) {
		T value = T();
		for (auto& op : options) {
			if (op.name == name) {
				value = op.get_value<T>();
				break;
			}
		}
		return value;  
	}
	
	/*! Sets the Option with the given name to the given value.
		If the Option already exists in the list of options for this command,
		it changes the value of that option; Options are not strongly typed,
		as all values are piped into stringstreams. If the Option does not
		exist, it creates a new Option and adds it to the list of options.
		@tparam T - the typing of the value for this Option.
		@param name - the name of the Option.
		@param val - the value of the Option.
	 */
	template <typename T>
	void set_option(std::string name, T val) {
		bool found = false;
		for (auto& op : options) {
			if (op.name == name) {
				op.set_value<T>(val);
				found = true;
				break;
			}
		}
		if (!found) {
			std::ostringstream os;
			os << val; 
			options.push_back(Option(name, os.str()));
		} 
	} 
	
	/*! Returns true if the Command has an Option with the given name, and false otherwise.
		@param name - the name of the Option.
		@return Boolean corresponding to whether the Option was found.
	 */
	bool is_option_set(std::string name) {
		bool found = false;
		for (auto& op :options) {
			if (op.name == name) {
				found = true;
				break;
			}
		}
		return found; 
	}
	
	int molecule_id; ///< The ID of the Molecule this command was called on.
	std::string name;  ///< The name of the Command.
};

/*! \struct Construct 
	\brief An abstraction of any input contained in curly brackets. 

	Any input of the form 
		name = { [CONTENT] }
	will be read into a construct. This can itself contain nested subconstructs,
	or just lines of text. These will then be parsed as either Constructs, Commands, 
	or Options, in a recursive descent algorithm.
 */
struct Construct {
	std::vector<Construct> subconstructs; ///< Collection of nested Constructs.
	std::vector<std::string> content;  ///< Input within the Construct that is not itself a Construct.
	std::string name; /// The name of the Construct.
	
	/// Creates an empty Construct with name "default".
	Construct() : name("default") { } 
	
	/*! Creates a Construct with the given name by parsing the text into subconstructs and content. 
		@param lines - the lines of strings contained within the curly brackets.
		@param name - the name of the Construct.
	 */
	Construct(std::vector<std::string>& lines, std::string& name); 
	
	/// Copy constructor
	Construct(const Construct& other);
	
	/*! Recursively parses the given text into Constructs and other input.
		@param lines - the text to be parsed.
	 */
	virtual void parse(std::vector<std::string>& lines);
};

/*! \class ProgramController
	\brief A container that reads, stores, and executes the input to the program.

	This takes the input and recursively parses it into Constructs, Commands, and
	Options. It currently assumes that base Constructs (i.e. those not contained
	within another Construct) are Molecules, and executes all commands following
	a Molecule on that Molecule until a second Molecule is reached. 

	The ProgramController is designed to be unique for each input file, and passed
	around the Program as a SharedPC (a shared pointer of a single ProgramController).

	It also contains the Logger, which handles all output.
 */
class ProgramController : public std::enable_shared_from_this<ProgramController> {
private:
	std::vector<Option> global_options; ///< Base Options for the whole program, i.e. those not contained in a Command or Construct.
	std::vector<Command> commands; ///< All Commands found in the input, labelled by Molecule index (position in constructs).
	std::vector<Construct> constructs;  ///< All base Constructs in the input, assumed to be Molecules.
	std::map<std::string, std::function<void(Command&, SharedMolecule)>> command_list; ///< A map pointing each command name to a calling function.
	
	std::shared_ptr<Fock> focker; ///< Pointer to the Fock object for the current Molecule/Command.
	std::shared_ptr<SCF> hf; ///< Pointer to the SCF object for the current Molecule/Command.
	std::shared_ptr<MP2> mp2obj;  ///< Pointer to the MP2 object for the current Molecule/Command.
	std::shared_ptr<IntegralEngine> ints; ///< Pointer to the IntegralEngine for the current Molecule.
	
	std::string jobname; ///< Name of the job, corresponding to the input file name with the extension removed.
	
	bool done_hf; ///< Flag for if an HF calculation has been done - further methods cannot be done otherwise.
	bool done_transform; ///< Flag for if an integral transformation has been done, necessary for CC methods.
	
	/*! A utility function for all methods that require an MP2-like initialisation (e.g. integral transformation)/
	 	@param mp2obj - a reference to the MP2 object to perform the integral transformation etc. on.
		@param hf - a reference to the SCF object containing the MO coefficients.
		@param calc - true if the MP2 energy is to be calculated directly, false otherwise.
	 */
	void runmp2(MP2& mp2obj, SCF& hf, bool calc); 
	
public:
	
	Logger log; ///< The Logger for the program, handling all output.
	
	/*! Creates a new ProgramController given streams for the input and output files, and the name of the program.
		@param input - the input file stream.
		@param output - the main output file stream.
		@param err - the error output stream.
		@param name - the jobname.
	 */
	ProgramController(std::ifstream& input, std::ofstream& output, std::ostream& err, std::string& name); 
	
	/// Overloaded copy operator.
	ProgramController& operator=(const ProgramController& other);
	
	/// Returns the jobname.
	std::string& getName() { return jobname; }
	
	/*! Returns the value of a global program Option, if found, as the specified type.
		@tparam T - the return type of the value for the desired Option.
		@param name - the name of the Option.
		@return The value of the global Option with the given name if it can be found. Otherwise, the default value for type T.
	 */
	template <typename T>
	T get_option(std::string name) {
		T value = T(); 
		for (auto& op : global_options) {
			if (op.name == name) {
				value = op.get_value<T>();
				break;
			}
		}
		return value;
	}
	
	/*! Sets the global Option with the given name to the given value.
		It first searches to see if the Option has already been defined and changes the value;
		if not a new global Option is created with the given name and value.
		@param name - the name of the Option.
		@param val - the value of the Option.
	 */
	template <typename T>
	void set_option(std::string name, T val) {
		bool found = false;
		for (auto& op : global_options) {
			if (op.name == name) {
				op.set_value<T>(val); 
				found = true;
				break;
			}
		}
		if (!found) {
			std::ostringstream os;
			os << val; 
			global_options.push_back(Option(name, os.str()));
		}
	}
	
	/*! Returns true if a global Option with the given name exists.
	 	@param name - the name of the desired Option.
		@return True if a global Option with the given name exists, false otherwise.
	 */
	bool is_option_set(std::string name) {
		bool found = false;
		for (auto& op : global_options) {
			if (op.name == name) {
				found = true;
				break;
			}
		}
		return found; 
	}
	 
	/*! Recursively parses the input file into Constructs, Commands, and Options.
		@param input - the input file stream.
	 */
	void parse(std::ifstream& input); 
	
	/// Runs the program; need only be called once.
	void run(); 
	
	/*! Utility function that removes comments and whitespace from a line of input.
	 	@param line - reference to a line of input to be cleaned.
	 */  
	void cleanLine(std::string& line, bool to_lower = true); 
	
	void call_hf(Command& c, SharedMolecule m); ///< Calls RHF or UHF depending on spin multiplicity.
	void call_rhf(Command& c, SharedMolecule m); ///< Calls RHF.
	void call_uhf(Command& c, SharedMolecule m); ///< Calls UHF.
	void call_mp2(Command& c, SharedMolecule m); ///< Calls MP2, if an HF calculation has been done.
	void call_dfmp2(Command& c, SharedMolecule m); ///< Calls density-fitted MP2, if an HF calculation has been done.
	void call_rpa(Command& c, SharedMolecule m); ///< Calls RPA, if an HF calculation has been done.
	void call_ccsd(Command& c, SharedMolecule m); ///< Calls CCSD, if an HF calculation and integral transformation has been done.
	void call_ralmo(Command& c, SharedMolecule m); ///< Calls restricted ALMO SCF.
	void call_ualmo(Command& c, SharedMolecule m); ///< Calls unrestricted ALMO SCF.
	void call_optg(Command& c, SharedMolecule m); ///< Calls a geometry optimization (currently only RHF or MP2)  
	void call_optx(Command& c, SharedMolecule m); ///< Calls a basis function exponent optimization (currently onl RHF or MP2)
	void call_nuctest(Command& c, SharedMolecule m); ///< Do not use - under development.
	
};

using SharedPC = std::shared_ptr<ProgramController>; ///< A shared pointer to a ProgramController.

#endif
