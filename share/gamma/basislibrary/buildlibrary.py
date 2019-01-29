import parsebasis as pb
import makelist as ml 

def build_library(name, parser=pb.parse, writer=pb.write_basis):
    file_list, name_list = ml.write_list_file(name)
    prefix = "./bases/build/" + name
    
    for i in range(len(file_list)):
        filename = file_list[i]
        basisname = [n for n, v in name_list.items() if v == i][0]
        
        with open(filename, 'r') as f:
            atoms = parser(f)
            writer(atoms, prefix, basisname, i)

build_library("basis")
build_library("jkfit")
build_library("mp2fit")
build_library("ecp", parser=pb.parse_ecp, writer=pb.write_ecp_basis)
    
    
