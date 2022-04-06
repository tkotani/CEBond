#!/usr/bin/env python3
import os
import re
import sys
#import pickle
import warnings
import subprocess
import numpy as np
warnings.filterwarnings('ignore')

import logging 
from pymatgen.io.cif import CifParser
from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import MultiWeightsChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import  SimplestChemenvStrategy, MultiWeightsChemenvStrategy

def adjacent_table1(cif_file):
        ### Setup the local geometry finder ###
        checktag = True
        lgf = LocalGeometryFinder()
        lgf.setup_parameters(centering_type='centroid', include_central_site_in_centroid=True)
        try:
                parser = CifParser(cif_file)
                try:
                        struct = parser.get_structures()[0]
                except:
                        print('error parser.get_structures()[0]')
                        pass
                        checktag = False
        except:
                pass
                
        try:
                lgf.setup_structure(structure=struct)
        except:
                pass
                print('error lgf.setup_structure(structure)')
                #checktag = ' error lgf.setup_structure(structure) '
                checktag = False
        
        return checktag


def adjacent_table2(cif_file,cifnumname,outputdir):
        ### Setup the local geometry finder ###
        lgf = LocalGeometryFinder()
        lgf.setup_parameters(centering_type='centroid', include_central_site_in_centroid=True)
        
        cif_text=subprocess.getoutput('cat {cif}'.format(cif=cif_file))
        try:
                space_group_IT_number = cif_text.split('_space_group_IT_number')[1].split('\n')[0].strip()
        except:
                space_group_IT_number = None
        
        parser = CifParser(cif_file)
        struct = parser.get_structures()[0]
        lgf.setup_structure(structure=struct)
        
        ### Get the StructureEnvironments ###
        try:    
                se = lgf.compute_structure_environments(maximum_distance_factor=1.41,only_cations=True)
        except:
                print('error lgf.compute_structure_environments(maximum_distance_factor=1.41,only_cations=True)')
                return 0
        
        ### Get lightstructure emvironment ###
        try:
                ###  strategy  ###   
                strategy = MultiWeightsChemenvStrategy.stats_article_weights_parameters()
                lse = LightStructureEnvironments.from_structure_environments(strategy=strategy,structure_environments=se)
        except:
                print('error LightStructureEnvironments.from_structure_environments(strategy=strategy,structure_environments=se)')
                return 0
        
        ### print ###
        with open( outputdir + "/" + "ntable_{0}.dat".format(cifnumname),"w") as file2:
                print('cif name=',cif_file,file=file2)
                print(cif_file,file=file2)
                print(lse.structure,file=file2)
                print(file=file2)
                print('space_group_IT_number',space_group_IT_number,file=file2)
                print('\n*** lse ***\n',file=file2)
                for isite in range(len(lse.structure)):
                        print("\n**isite=",isite,file=file2)
                        
                        #print(se.structure[isite].species_and_occu,":",lse.structure.sites[isite].coords,file=file2)
                        print(se.structure[isite].species,":",lse.structure.sites[isite].coords,file=file2)
                        
                        isite_coords=lse.structure.sites[isite].coords
                        print(file=file2)
                        if lse.coordination_environments[isite] is None:
                                continue 
                        for i in range(len(lse.coordination_environments[isite])):
                                print('ce_symbol  =',lse.coordination_environments[isite][i]['ce_symbol'],file=file2)
                                print('ce_fraction= ',lse.coordination_environments[isite][i]['ce_fraction'],file=file2)
                                cn=len(lse.neighbors_sets[isite][i].all_nbs_sites_indices)
                                _distance=[]
                                print("species:coords:index",file=file2)
                                for inbs in range(len(lse.neighbors_sets[isite][i].neighb_sites_and_indices)):      #neib loop                   
                                        coords=lse.neighbors_sets[isite][i].neighb_sites_and_indices[inbs]['site'].coords
                                        index=lse.neighbors_sets[isite][i].neighb_sites_and_indices[inbs]['index']
                                        coords_checkindex=lse.structure.sites[index].coords

                                        _distance.append(np.linalg.norm(isite_coords - coords))
                                        species=lse.structure.sites[index]._species
                                        print(species,":",coords,":",index,file=file2)
                                
                                if not set(se.neighbors_sets[isite][cn][0].distances)==set(_distance):
                                        print("!!!atention to nb coordination!!!",file=file2)
                                        print("        ",se.neighbors_sets[isite][cn][0].distances,file=file2)
                                        print("        ",_distance,file=file2)
                                print(file=file2)
                print('END\n',file=file2)

        return 0

def main():
        cwd = os.getcwd()
        anion = sys.argv
        ndir= cwd+'/neighbour_'+anion[1]
        outputdir=ndir
        if not os.path.isdir(ndir): os.mkdir(ndir)
        print('OutputDir=',outputdir)
        ciflist = open("list_{0}.dat".format(anion[1]),"r").read().split('\n')
        os.chdir('./COD/'+anion[1])
        cwd = os.getcwd()
        print(cwd)
        for i in ciflist: 
                try:
                        cif_file = i.split()[0]
                except:
                        break
                print(cif_file,' ------------------------------ ')
                cifnumname = re.findall('\/(\d+).cif',cif_file)[0]
                if adjacent_table1(cif_file):
                        adjacent_table2(cif_file,cifnumname,outputdir)


if __name__ == '__main__':
        main()

