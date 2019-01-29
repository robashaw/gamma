from lxml import etree

angular_shells = ["s", "p", "d", "f", "g", "h", "i"]

class Shell:
    def __init__(self, lval = 0, nexp = 0, nbfs = 0):
        self.lval = lval
        self.powers = []
        self.nexp = nexp
        self.exps = []
        self.contr = []

class Atom:
    def __init__(self, name="X", nshells = 0):
        self.name = name
        self.ncore = 0
        self.maxl = 0
        self.nshells = nshells
        self.shells = []
        
def parse(file):
    atoms = []
    atomnames = {}
    lines = file.readlines()
    nlines = len(lines)
    linenumber = 0
    atomnumber = 0
    
    while linenumber < nlines:
        line = lines[linenumber].strip()
        tokens = line.split(',')
        for token in tokens: 
            token = token.replace(' ', '')
            
        if (len(tokens) > 1):
            shell = tokens[0].replace(' ', '')
             
            if (shell in angular_shells):
                atom_name = tokens[1].replace(' ', '').lower()
                current_atomnumber = atomnumber
                if atom_name in atomnames:
                    current_atomnumber = atomnames[atom_name]
                else:
                    atomnames[atom_name] = atomnumber
                    new_atom = Atom(name=atom_name)
                    atoms.append(new_atom)
                    atomnumber += 1
                current_atom = atoms[current_atomnumber]
                
                exps = tokens[2:]
                nexp = len(exps)
                
                end_of_shell = False
                nfuncs = 0
                
                while not end_of_shell:
                    linenumber += 1
                    line = lines[linenumber].strip()
                    tokens = line.split(',')
                    for token in tokens: 
                        token = token.replace(' ', '')
                        
                    if len(tokens) > 1:
                        if tokens[0].replace(' ', '') is "c":
                            xix = tokens[1].replace(' ', '').split('.')
                            start = int(xix[0])-1
                            end = int(xix[1])-1

                            new_shell = Shell(lval=shell, nexp=end-start+1, nbfs=1)
                            new_shell.exps = exps[start:end+1]
                            new_shell.contr = tokens[2:]
                            
                            current_atom.shells.append(new_shell)
                            current_atom.nshells += 1
                        else:
                            end_of_shell = True
                            linenumber -= 1
                    else:
                        end_of_shell = True
                        linenumber -=1           
        linenumber += 1
    
    return atoms
                
def write_basis(atoms, prefix, name, index):
    filename = prefix + str(index) + ".xml"
    with open(filename, 'wb') as new_file:
        root = etree.Element("root", name=name)
        tree = etree.ElementTree(root)
        
        for atom in atoms:
            child = etree.SubElement(root, atom.name, nshells = str(atom.nshells))
            
            for shell in atom.shells:
                schild = etree.SubElement(child, "Shell", lval=shell.lval, nexp=str(shell.nexp))
                
                for i in range(shell.nexp):
                    try:
                        xchild = etree.SubElement(schild, "xc", x=shell.exps[i], c=shell.contr[i])
                    except:
                        print("ERROR in " + filename + " (" + name + ") : atom " + atom.name + ", shell type " + shell.lval)
                        print("Expected no. of exps/coeffs:", shell.nexp)
                        print("Actual no. of exps: ", len(shell.exps))
                        print("Actual no. of coeffs: ",  len(shell.contr))
        tree.write(new_file, pretty_print = True)

def parse_ecp(file):
    atoms = []
    atomnames = {}
    lines = file.readlines()
    nlines = len(lines)
    linenumber = 0
    atomnumber = 0
    
    while linenumber < nlines:
        line = lines[linenumber].strip()
        tokens = line.split(',')
        for token in tokens: 
            token = token.replace(' ', '')
            
        if (len(tokens) > 2):
            if (tokens[0].lower() == "ecp"):
                atom_name = tokens[1].lower()
                new_atom = Atom(name=atom_name)
                
                new_atom.ncore = int(tokens[2])
                new_atom.maxl = int(tokens[3])
                
                for i in range(new_atom.maxl+1):
                    linenumber += 1
                    line = lines[linenumber].strip()
                    
                    tokens = line.split(';')
                    for token in tokens:
                        token = token.replace(' ', '')
                        
                    if(len(tokens) > 1):
                        nprims = int(tokens[0])
                        l = i-1
                        if (i==0):
                            l = new_atom.maxl
                        new_shell = Shell(lval=l, nexp=nprims, nbfs=1)
                        
                        for token in tokens[1:]:
                            subtokens = token.split(',')
                            if (len(subtokens) == 3):
                                new_shell.powers.append(subtokens[0])
                                new_shell.exps.append(subtokens[1])
                                new_shell.contr.append(subtokens[2])
                        
                        new_atom.shells.append(new_shell)
                
                atoms.append(new_atom)
                
        linenumber += 1
    
    return atoms

def write_ecp_basis(atoms, prefix, name, index):
    filename = prefix + str(index) + ".xml"
    with open(filename, 'wb') as new_file:
        root = etree.Element("root", name=name)
        tree = etree.ElementTree(root)
        
        for atom in atoms:
            child = etree.SubElement(root, atom.name, ncore = str(atom.ncore), maxl=str(atom.maxl))
            
            for shell in atom.shells:
                schild = etree.SubElement(child, "Shell", lval=str(shell.lval), nexp=str(shell.nexp))
                
                for i in range(shell.nexp):
                    try:
                        xchild = etree.SubElement(schild, "nxc", n=shell.powers[i], x=shell.exps[i], c=shell.contr[i])
                    except:
                        print("ERROR in " + filename + " (" + name + ") : atom " + atom.name + ", shell type " + shell.lval)
                        print("Expected no. of powers/exps/coeffs:", shell.nexp)
                        print("Actual no. of powers: ", len(shell.powers))
                        print("Actual no. of exps: ", len(shell.exps))
                        print("Actual no. of coeffs: ",  len(shell.contr))
        tree.write(new_file, pretty_print = True)
    
