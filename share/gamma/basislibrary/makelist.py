import os
from lxml import etree

def make_file_list(extension=".basis"):
    file_list = []
    for root, dirs, files in os.walk("./bases"):
        for file in files:
            if file.endswith(extension):
                filepath = os.path.join(root, file)
                file_list.append(filepath)
    return file_list

def make_name_list(file_list):
    name_list = {}
    file_ctr = 0
    for file in file_list:
        with open(file, 'r') as f:
            first_line = f.readline().strip()
            names = first_line.split(',')
            
            for name in names:
                clean_name = name.replace(' ', '')
                name_list[clean_name] = file_ctr
        file_ctr += 1
    return name_list

def write_list_file(name):
    extension = "." + name
    file_list = make_file_list(extension)
    name_list = make_name_list(file_list)
    
    list_file = name + ".list"
    with open(list_file, 'wb') as lf:
        root = etree.Element("root", name=name)
        tree = etree.ElementTree(root)
        
        for key in name_list:
            child = etree.SubElement(root, "Item", key=key, value=str(name_list[key]))
        
        tree.write(lf, pretty_print=True)
    return file_list, name_list   
