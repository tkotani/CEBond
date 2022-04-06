#!/usr/bin/env python3
"""
usage:
>CifFilter.py O
"""
import os
import re
import sys
import pickle
import subprocess

import logging 
import pymatgen as mg

import skiplist as sk
skip_list= sk.skiplist

def check_search_anion(ciffilelist): #,search_anion,eliminate_anion_list):
        total_ciffile = []
        for i in ciffilelist:
                cifname=i.rsplit("/",1)[1].split(".")[0]
                if cifname in skip_list:
                        print("        skip in skiplist.py ",cifname+'.cif')
                        continue
                try:
                        structure = mg.Structure.from_file(i)
                        composition = structure.species_and_occu
                        #print(composition)
                except:
                        print("        error pymatgen.Structure ",cifname+'.cif')
                        #error_structure_ciffile.append(i)
                        continue
                
                skip=False
                for j in composition:
                        occu = float(re.split(r'([a-zA-Z+-]+)',str(j))[-1]) #get last number 
                        #print(j,occu)
                        if( occu<0.9 and occu>0.1):
                                #print(composition)
                                print("        skip non-stoichiometic ",cifname+'.cif')
                                skip=True
                                break
                if skip: continue
                print(i)
                total_ciffile.append("{} : {}".format(str(i),str(structure.formula)))
        return total_ciffile

def main():
        cwd = os.getcwd()   #Working dirctory
        cod_path = './COD/'+sys.argv[1]   #Cod dirctory
        search_anion = sys.argv[1]   #Target anion
        
#eliminate_anion_list = ['O','N','F','P','S','Cl','As','Se','Br','Sb','Te','I','Bi','Po','At','H','C']
#eliminate_anion_list.remove(search_anion) #Eliminate_anion_list : anion list how except target anion 
        #print(eliminate_anion_list)
        
        os.chdir(cod_path)
        ciffilelist_ = subprocess.getoutput("find . -name '*.cif'|sort")  #search cif files
        ciffilelist = ciffilelist_.split('\n')  
        TargetCifFile = check_search_anion(ciffilelist) #,search_anion,eliminate_anion_list)
        os.chdir(cwd)
        file1name = "list_{0}.dat".format(search_anion)
        with open(file1name,"w") as file1:
                for i in TargetCifFile: #[0]:
                        print(i,file=file1)
                
        #with open("TargetCifFileList_{0}.pickle".format(search_anion),"wb") as fwb:
        #           pickle.dump(TargetCifFile[0],fwb)

if __name__ == '__main__':
        main()

