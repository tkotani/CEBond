#!/usr/bin/env python3
# Ver 3.6 2021 03 change layout
# Ver 3.3 2020 04 13 new color
# Ver 3.2            check High press Mg cig 
#                       and adapt hatch color(test)
#                       and output analized list
# Ver 3.1 2020 02 19 adapt file average analize option
# Ver 3.0 2020 02 new press and temp filter 
# Ver2.5 2019 12 16 correct bond weight
# Ver2.4 2019 12 11 adapts log scale option 
# distance histgram ce color,shannon arrow ,
# valence analize range R site 

import warnings
warnings.filterwarnings('ignore')
import sys,os,subprocess
import numpy as np
import re

import shutil

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import module_valence_bond_analizer_RBVS as vba
import graf_module as gm

###  "keyword=aaa" take aaa
def keycut(content,keyword):   
         return content.split(keyword)[1].split('\n')[0]

def make_bondhist(figset_name,anion,average_analize): 
        valence_cut_point=[]

        ###### no valence sets ################
        ### alcari set
        if figset_name=="alkali":
                alcari=True
                cation_list=["Li","Na","K","Rb","Cs"]
                valence_cut_point=[]
                fit_valence=[1]
                #xlim=[1.6,3.8]  #for alcari or alcari arth
                xlim=[1.5,3.7]  #for alcari or alcari arth

        ### alcari erth set
        if figset_name=="alkali_arth":
                alcari=True
                cation_list=["Be","Mg","Ca","Sr","Ba"]
                valence_cut_point=[]
                fit_valence=[2]
                #xlim=[1.4,3.6]  #for alcari arth
                xlim=[1.5,3.7]  #for alcari arth

        #cation_list=["Mg"]
        #xlim=[1.0,2.5]
        high_press_check=""
        analized_list=""

        ### transitional metal set  (No valence)
        if figset_name=="3d_1":
                alcari=False
                cation_list=["Sc","Ti","V","Cr","Mn"]
                fit_valence=[]
                xlim=[1.4,2.5] #for transitonal metal

        if figset_name=="3d_2":
                alcari=False
                cation_list=["Fe","Co","Ni","Cu","Zn"]
                fit_valence=[]
                xlim=[1.4,2.5] #for transitonal metal

        if figset_name=="4d_1":
                alcari=False
                cation_list=["Y","Zr","Nb","Mo","Tc"]
                valence_cut_point=[]
                fit_valence=[]
                xlim=[1.4,2.5] #for transitonal metal
        if figset_name=="4d_2":
                alcari=False
                cation_list=["Ru","Rh","Pd","Ag","Cd"]
                valence_cut_point=[]
                fit_valence=[]
                xlim=[1.4,2.5] #for transitonal metal
        
        if figset_name=="abst":
                alcari=False
                cation_list=["V","Cr","Mn"]
                valence_cut_point=[]
                fit_valence=[]
                xlim=[1.4,2.4] #for transitonal metal

        if figset_name=="valence":
                alcari=False
                
                if os.path.exists('./contena/VBA_dat'):
                        shutil.rmtree('./contena/VBA_dat')
                
                cation_type = sys.argv[2]
                if cation_type == "V":
                        cation_list=["V"]
                        #valence_cut_point=[2.5,4]
                        fit_valence=[2,3,4,5]
               
                if cation_type == "Cu":
                        cation_list=["Cu"]
                        fit_valence=[1,2]
                
                if cation_type == "Fe":
                        cation_list=["Fe"]
                        fit_valence=[2,3,6]
                           
                #xlim=[1.4,2.4]
                xlim=[1.4,2.5]


        ###### split for valence sets ##############
        """ old set
        #alcari=False

        #cation="Li"
        #fit_valence=[1]
        #valence_cut_point=[]

        #cation="Na"
        #fit_valence=[1]
        #valence_cut_point=[]


        #cation="Cu"
        #fit_valence=[1,2]
        #valence_cut_point=[1.5]

        #cation="V"
        #fit_valence=[2,3,4]
        #valence_cut_point=[2.5,3.5]

        #cation="Ni"
        #fit_valence=[2,3]
        #valence_cut_point=[2.5]

        #cation="Co"
        #fit_valence=[2,3,4]
        #valence_cut_point=[2.5,3.5]

        #cation="Zn"
        #fit_valence=[2]
        #valence_cut_point=[]


        ###data="neib_table_8.0_Fe"
        #cation_list=["Fe"]
        #fit_valence=[2.0,2.2,2.4,2.6,2.8,3.0]
        #valence_cut_point=[2.1,2.3,2.5,2.7,2.9]
        #fit_valence=[2,3]
        #valence_cut_point=[2.5]
        #valence_cut_point=[]   #[]  list is None  do not split for valence
        """

        #### !!!check alcari=False  ######
        #alcari=False
        #xlim=[1.4,2.5] #for transitonal metal

        """3d metal set"""
        #cation_list=["Sc"]
        #fit_valence=[3]

        #cation_list=["Ti"]
        #fit_valence=[2,3,4]

        #cation_list=["V"]
        #fit_valence=[2,3,4,5]
        #fit_valence=[round(2.8+i*0.2,1) for i in range(2+5*2)]

        #cation_list=["Cr"]
        #fit_valence=[2,3,4,6]

        #cation_list=["Mn"]
        #fit_valence=[2,3,4,6,7]

        #cation_list=["Fe"]
        #fit_valence=[2,3,6]
        #fit_valence=[1,2,3,4,5,6]
        #fit_valence=[round(1.8+i*0.2,1) for i in range(7)]

        #cation_list=["Co"]
        #fit_valence=[2,3]

        #cation_list=["Ni"]
        #fit_valence=[2,3]

        #cation_list=["Cu"]
        #fit_valence=[1,2]
        #valence_cut_point=[1.7]

        #cation_list=["Zn"]
        #fit_valence=[2]

        """4d metal set"""
        #cation_list=["Y"]
        #fit_valence=[3]

        #cation_list=["Zr"]
        #fit_valence=[4]

        #cation_list=["Nb"]
        #fit_valence=[3,5]

        #cation_list=["Mo"]
        #fit_valence=[1,2,3,4,5,6]

        #cation_list=["Tc"]
        #fit_valence=[4,5,6,7]

        #cation_list=["Ru"]
        #fit_valence=[2,3,4,5,6,7]

        #cation_list=["Rh"]
        #fit_valence=[3]

        #cation_list=["Pd"]
        #fit_valence=[2,4]

        #cation_list=["Ag"]
        #fit_valence=[1]

        #cation_list=["Cd"]
        #fit_valence=[2]

        ### fit_valence=[2,3]  >>>   valence_cut_point=[2.5]
        if len(fit_valence)>=2 and not valence_cut_point:
                valence_cut_point=[(fit_valence[v]+fit_valence[v+1])/2 for v in range(len(fit_valence)-1)]

        ######################### <Options> ######################## 
        nb_only_O=True  #take only cation-O distance
        R=5             #range for valence Analize
        class_width=0.01
        #class_width=0.025
        limit_site=500  #None  #None or Number
        cn7_base=True
        log_scale=False #log scale option
        average_analize=False#True
        value_log=False
        valencegraf=True
        ######################### anion sets #######################
        if anion=="O":
                R_O=1.4         #eR for O  needed shannon arrow
                #anion='O'
                list_file="list_O_P-T"
        if anion=="N":
                #anion='N'
                #fit_valence=[]
                R_O=1.46
                list_file="list_N_P-T"
        if anion=="F":
                #anion='F'
                #fit_valence=[]
                R_O=1.33 ## (R_F)
                list_file="list_F_P-T"

        #import text saved chem data
        #with open(data) as l:
        #        contents=l.read()

        ### for graf
        contena={}
        max_distance=10#5
        class_width  =0.01                                                                      
        loop_count=0
        loop_count1=0
        ### main coords ###
        #for one_cif_data in contents.split("cif name= ")[1:]:

        ##Ver 3.0 input dat file
        import module_text_contena as mt
        inport_list=mt.text_input(list_file).split("\n")
        dir_list=[item.split(" ")[1] for item in inport_list if not item==""]

        if valencegraf==True:
                vmax_distance=10
                vclass_width=0.1
                vcontena=[0]*int(vmax_distance/vclass_width)

        ### input skip list
        skip_list1="./contena/errors/high_press_Mg_1.40<distance<1.42"
        skip_list=mt.text_input(skip_list1)

        for cation in cation_list:
                ### cif loop
                for cif_path in dir_list:
                        print("\n\n")
                        print('***************** cif ***************************')

                        # input nb table
                        cif_num=cif_path.strip(".cif").split("/")[1]

                        if cif_num in skip_list:
                                print("skip",cif_num)
                                continue

                        try:one_cif_data=mt.text_input("./contena/neib_tables/ntable_"+cif_num+".dat")              #input nbtable(dat) 
                        except:
                                print("!!!none dat file")
                                continue

                        # count
                        if one_cif_data=="\n":
                                loop_count1+=1
                                continue

                        loop_count+=1
                        cif_weight=0
                        _cif_weight=0
                        high_press_switch=False

                        struct_data=one_cif_data.split("*** lse ***")[0]
                        lse_data   =one_cif_data.split("*** lse ***")[1]
                        cifname=keycut(one_cif_data,"cif name=")
                        print(cifname)
                        cifname=re.search(r'\d{7}',cifname).group()
                        
                        valence_list=[]
                        center_symbols=[]
                        center_points=[]

                        #in advance make center data list       >>>center_symbols , center_points
                        for site_data in lse_data.split('\n\n\n'):  
                                if site_data=='': continue  
                                center_data=site_data.split("\n\n")[0].split("\n")[1]
                                
                                if center_data=='' :continue

                                center_symbol=center_data.split(":")[0]
                                center_point =center_data.split(":")[1]
                                
                                center_point =[float(point) for point in center_point.strip().strip("[").strip("]").split()]

                                center_symbols.append(center_symbol)
                                center_points.append(center_point)



                        #### !!!! test limit a number of site
                        if limit_site:
                                if len(center_symbols)>limit_site:
                                        print("!!! nam site>{l}".format(l=limit_site))
                                        continue

                        ### check fractional oqupancy site
                        #print("<check fractional oqupancy site>")
                        have_fractional_site=None
                        for site in center_symbols:
                                atoms=re.findall(r'[A-Z][a-z]*',site)
                                #print("atoms",atoms)
                                #print(len(atoms))
                                if len(atoms)>1.1:
                                       have_fractional_site=True
                        if have_fractional_site :
                                print("!!!! fractional oqupancy site  !!!!")
                                continue  ## skip this cif by have a fractional occupancy site   
                  
                        center_symbols=[re.match(r'[A-Z][a-z]*',atom).group() for atom in center_symbols]  ## ex) Fe2+1 >> Fe

                        n_cation=0  # a number of analized cation 
                        for atom in center_symbols:
                                if re.match(r'[A-Z][a-z]*',atom).group() == cation:
                                        n_cation+=1
                        print('n_cation= ',n_cation)

                        # search cation
                        if n_cation==0:
                                print("! no cation",cation )
                                continue       

                        ### valence block ###
                        if valence_cut_point:
                                ##### take valence list    method 1
                                print(struct_data)
                                abc=[float(d) for d in keycut(struct_data,"abc   :").split()]
                                angles=[float(d) for d in keycut(struct_data,"angles:").split()]

                                #site_data=re.split(r'\s\s-+\n',struct_data)[1].split("\nvalence=")[0].split("\n")
                                site_data=re.split(r'\s\s-+\n',struct_data)[1].split("\nspace_group_IT_number")[0].split("\nvalence=")[0].split("\n")

                                site_list=[]
                                fcoord_list=[]
                                for sline in site_data:
                                        if sline=='':continue
                                        i=sline.split()
                                        site=re.match(r'[A-Z][a-z]*',i[1]).group()
                                        coord=[float(c) for c in i[2:5]]
                                        site_list.append(site)
                                        fcoord_list.append(coord)

                                try:   
                                        content = mt.text_input('./contena/VBA_dat/{icif}_out.dat'.format(icif=cifname))  ## cash out.dat exist
                                        content=None  #reset 
                                        print("chash exist  skip trance analize")
                                        valence_list=vba.BVS_Analizer_R_2(cifname,R,site_list)
                                except:
                                        #trance lattice vectol
                                        print("chash none")
                                        print("start lattice vector")
                                        a,b,c=vba.trance_coordination_f_to_cartesian(abc[0],abc[1],abc[2],angles[0],angles[1],angles[2])
                                        print("finish lattice vector")

                                        #trance fcoods to cartegian coords
                                        print("start trance coords")
                                        coord_list=[]
                                        for fcoord in fcoord_list:
                                             coord=[None]*3
                                             coord[0]=a[0]*fcoord[0] +b[0]*fcoord[1] +c[0]*fcoord[2]
                                             coord[1]=a[1]*fcoord[0] +b[1]*fcoord[1] +c[1]*fcoord[2]
                                             coord[2]=a[2]*fcoord[0] +b[2]*fcoord[1] +c[2]*fcoord[2]

                                             coord_list.append(coord)
                                        
                                        print("latice =",a,b,c)
                                        print("finish trance coords")

                                        valence_list=vba.BVS_Analizer_R_2(cifname,R,site_list,a,b,c,coord_list)     
                                print("\nvalence")
                                for i in range(len(site_list)):
                                        print(site_list[i],valence_list[i])     #print site and valence
                        else: valence_list=None
                        
                        ##### count cations with O site in neighbor 
                        invalid_n_cations=0     #a num of no anion bond cation
                        for site_data in lse_data.split('\n\n\n'):  #site data loop
                                if site_data=='': continue

                                center_data=site_data.split("\n\n")[0].split("\n")[1]
                                if center_data=='' :continue

                                center_data=re.match(r'[A-Z][a-z]*',center_data).group()
                                if not center_data==cation:
                                        continue

                                #print("site_data\n",site_data)
                                #neigbor data each ce
                                #print(site_data.split("\n\n"))
                                #print("\n")  
              
                                if len(site_data.split("\n\n"))==1:
                                        invalid_n_cations+=1

                                for neib_data in site_data.split("\n\n")[1:]:
                                        if (neib_data=="") or (neib_data=='END'):continue
                                        #print("neib ",neib_data)
                                        #if neib_data=='END':continue


                                        ce_symbol=keycut(neib_data,"ce_symbol  = ")
                                        ce_fraction=float(keycut(neib_data,"ce_fraction=  ")   )

                                        neib_point_data=neib_data.split("species:coords:index\n")[1]     

                                        neighbors=[]
                                        #neighbors_positions=[] 
                                        for one_neib_point in neib_point_data.split("\n"):
                                                if one_neib_point=='!!!atention to nb coordination!!!':break
                                                nb_symbol=one_neib_point.split(":")[0]
                                                #nb_point =one_neib_point.split(":")[1]
                                                #nb_point =one_neib_point.split(":")[1].strip().strip('[').strip(']').split()
                                                #nb_point_list=[float(point) for point in nb_point]

                                                neighbors.append(nb_symbol)
                                                #neighbors_positions.append(nb_point_list)
                                        
                                        neighbors=[re.match(r'[A-Z][a-z]*',atom).group() for atom in neighbors]
                                        if not anion in neighbors:
                                                print("! No O in neib list")

                                                invalid_n_cations+=ce_fraction
                        #print("invalid_n_cations ",invalid_n_cations)
                        n_cation=n_cation-invalid_n_cations
                        #print("n cation",n_cation)

                        ### ver2.2
                        ### in addvance make nb list
                        """
                        for site_data in lse_data.split('\n\n\n'):  
                                if site_data=='': continue
                                center_data=site_data.split("\n\n")[0].split("\n")[1]
                                if center_data=='' :continue
                                center_symbol=center_data.split(":")[0]
                                print("\n<",center_symbol," site analize>")

                                center_point =center_data.split(":")[1]
                                center_point =[float(point) for point in center_point.strip().strip("[").strip("]").split()]  #center point list
                                
                                #neigbor data each ce
                                for neib_data in site_data.split("\n\n")[1:]:
                                        if (neib_data=="") or (neib_data=='END'):continue

                                        ce_symbol=keycut(neib_data,"ce_symbol  = ")
                                        ce_fraction=float(keycut(neib_data,"ce_fraction=  ")   )

                                        neib_point_data=neib_data.split("species:coords:index\n")[1]     

                                        neighbors=[]
                                        neighbors_positions=[] 
                                        for one_neib_point in neib_point_data.split("\n"):
                                                if one_neib_point=='!!!atention to nb coordination!!!':break

                                                nb_symbol=one_neib_point.split(":")[0]
                                                nb_point =one_neib_point.split(":")[1]
                                                nb_point =one_neib_point.split(":")[1].strip().strip('[').strip(']').split()
                                                nb_point_list=[float(point) for point in nb_point]

                                                neighbors.append(nb_symbol)
                                                neighbors_positions.append(nb_point_list)

                                        print()
                                        print(ce_symbol,ce_fraction) 
                                        print("neighbors=",neighbors)
                                        print("neighbors position=",neighbors_positions)
                        """
                        ### fit valence
                        fit_valence_list=[]
                        if valence_cut_point:
                                for i in range(len(valence_list)):
                                        for cut in range(len(valence_cut_point)):
                                                if valence_list[i] < valence_cut_point[cut]:
                                                        fit_valence_list.append(fit_valence[cut])
                                                        break
                                        if valence_list[i] > valence_cut_point[len(valence_cut_point)-1]:   
                                                fit_valence_list.append(fit_valence[len(fit_valence)-1])  
                                print("fit valence list= ",fit_valence_list) 
                        elif fit_valence:fit_valence_list=[fit_valence[0]]*len(center_symbols) 
                        else:
                                #fit_valence_list=[""]*len(valence_list)                 #do not split for valence
                                fit_valence_list=[""]*len(center_symbols)                 #do not split for valence
                                                
                                print("fit valence list= ",fit_valence_list)   
                        for i in range(len(fit_valence_list)):
                                if fit_valence_list[i]==1:
                                        fit_valence_list[i]=''
                        

                        ### take distance and analize weight point 
                        isite=-1
                        #loop site
                        for site_data in lse_data.split('\n\n\n'):  
                                if site_data=='': continue
                                isite+=1
                                if not center_symbols[isite] == cation:
                                        continue
                                print("\n",center_symbols[isite])

                                #loop ce
                                for neib_data in site_data.split("\n\n")[1:]:
                                        if (neib_data=="") or (neib_data=='END'):continue

                                        ce_symbol   = keycut(neib_data,"ce_symbol  = ")
                                        ce_fraction = float(keycut(neib_data,"ce_fraction=  "))


                                        #loop neibors
                                        neib_point_data=neib_data.split("species:coords:index\n")[1]     
                                        neighbors=[]
                                        neighbors_positions=[] 
                                        for one_neib_point in neib_point_data.split("\n"):
                                                if one_neib_point=='!!!atention to nb coordination!!!':break

                                                nb_symbol=one_neib_point.split(":")[0]
                                                nb_point =one_neib_point.split(":")[1]
                                                nb_point =one_neib_point.split(":")[1].strip().strip('[').strip(']').split()
                                                nb_point_list=[float(point) for point in nb_point]

                                                neighbors.append(nb_symbol)                     #[neighbor1,neighbor2,,,,,]
                                                neighbors_positions.append(nb_point_list)       #[[x,y,z],[x,y,z],,,,,]

                                        neighbors=[re.match(r'[A-Z][a-z]*',symbol).group() for symbol in neighbors]

                                        print()
                                        print(ce_symbol,ce_fraction) 
                                        print("neighbors=",neighbors)
                                        print("neighbors position=",neighbors_positions)

                                        ### calcurate a bond point for this ce
                                        if nb_only_O: #take a only cation-O distance
                                                ### method 2
                                                if average_analize: ### average bonds distance  ,one site to one plot ,one file weight ==1
                                                #if True:
                                                        bond_distance_list=[]
                                                        for i_neigb in range(len(neighbors)):
                                                                if not neighbors[i_neigb]==anion: #not O bond skip
                                                                        print(isite,'-',i_neigb,"bond is no anion bond")
                                                                        continue
                         
                                                                bond_distance  = np.sqrt( (neighbors_positions[i_neigb][0]- center_points[isite][0])**2
                                                                                         +(neighbors_positions[i_neigb][1]- center_points[isite][1])**2
                                                                                         +(neighbors_positions[i_neigb][2]- center_points[isite][2])**2)
                                                                bond_distance_list.append(bond_distance)
                                                        if not bond_distance_list: continue

                                                        average_bond_distance=sum(bond_distance_list)/len(bond_distance_list)
                                                        bond_weight = ce_fraction/n_cation

                                                        if 1.40<bond_distance<1.42:     
                                                                high_press_switch=True        

                                                        # input to contena
                                                        #ex. fit_valence=[]
                                                        
                                                        #if not fit_valence:                                           #for tarancitonal metal no valence cut ex. [Sc Ti V Cr Mn]
                                                        if not fit_valence or alcari==True:
                                                                cation_key=cation                                     #contena={'Sc':{"O6":[d,d,d,d,d,d,,]}}
                                                        else:
                                                                cation_key=cation + str(fit_valence_list[isite])+'+'  #contena={"Sc3+":{"O6":[d,d,d,d,d,d,,]}}
                                                        
                                                        #cation_key=cation 


                                                        if not cation_key in contena.keys():
                                                                contena[cation_key]={}
                                                        if not ce_symbol in contena[cation_key].keys():
                                                                contena[cation_key][ce_symbol]=[0]*int(max_distance/class_width)

                                                        if not value_log:
                                                                contena[cation_key][ce_symbol][gm.distanse_divider(average_bond_distance,class_width)] += bond_weight          #contain weight                                                                   
                                                        else:
                                                                contena[cation_key][ce_symbol][gm.distanse_divider_log(average_bond_distance,class_width)] += bond_weight          #contain weight                                                                   
                                                        cif_weight+=bond_weight
                                                        #bond_weight1= bond_weight        

                                                ### method 1
                                                else:# one bond to one plot , one cif weight is 1
                                                        for i_neigb in range(len(neighbors)):
                                                                if not neighbors[i_neigb]==anion: #not O bond skip
                                                                        print(isite,'-',i_neigb,"bond is no anion bond")
                                                                        continue
                         
                                                                bond_distance  = np.sqrt( (neighbors_positions[i_neigb][0]- center_points[isite][0])**2
                                                                                         +(neighbors_positions[i_neigb][1]- center_points[isite][1])**2
                                                                                         +(neighbors_positions[i_neigb][2]- center_points[isite][2])**2)


                                                                bond_weight = ce_fraction/neighbors.count(anion)/n_cation

                                                                

                                                                # input to contena
                                                                #ex. fit_valence=[]
                                                                
                                                                #if not fit_valence:                                           #for tarancitonal metal no valence cut ex. [Sc Ti V Cr Mn]
                                                                if not fit_valence or alcari==True:
                                                                        cation_key=cation                                     #contena={'Sc':{"O6":[d,d,d,d,d,d,,]}}
                                                                else:
                                                                        cation_key=cation + str(fit_valence_list[isite])+'+'  #contena={"Sc3+":{"O6":[d,d,d,d,d,d,,]}}
                                                                
                                                                #cation_key=cation 

                                                                if not cation_key in contena.keys():
                                                                        contena[cation_key]={}
                                                                if not ce_symbol in contena[cation_key].keys():
                                                                        contena[cation_key][ce_symbol]=[0]*int(max_distance/class_width)

                                                                if not value_log:
                                                                        contena[cation_key][ce_symbol][gm.distanse_divider(bond_distance,class_width)] += bond_weight          #contain weight                                                                   
                                                                else:
                                                                        contena[cation_key][ce_symbol][gm.distanse_divider_log(bond_distance,class_width)] += bond_weight          #contain weight                                                                   
                                                                
                                                                cif_weight+=bond_weight
                                                                         
                                                                #bond_weight2= bond_weight 
                                                #if not round(bond_weight1,8)==round(bond_weight2,8):
                                                #        print("bond_weight=",bond_weight1,bond_weight2)
                                                #        #sys.exit() 

                                        else:   # take all distances
                                                for i_neigb in range(len(neighbors)):
                                                        bond_distance  = np.sqrt( (neighbors_positions[i_neigb][0]- center_points[isite][0])**2
                                                                                 +(neighbors_positions[i_neigb][1]- center_points[isite][1])**2
                                                                                 +(neighbors_positions[i_neigb][2]- center_points[isite][2])**2)
                                                        bond_weight = ce_fraction/len(neighbors)/n_cation
                                                        if not cation in contena.keys():
                                                                contena[cation]={}
                                                        if not ce_symbol in contena[cation].keys():
                                                                contena[cation][ce_symbol]=[0]*int(max_distance/class_width)
                                                        contena[cation][ce_symbol][gm.distanse_divider(bond_distance,class_width)] += bond_weight          #contain weight

                        #if not round(cif_weight,8)==round(_cif_weight,8):
                        #        print("cif_weight",cif_weight,_cif_weight) 
                        #        sys.exit()

                        if cif_weight==0 : 
                                print('no weight')
                                #sys.exit()
                                continue              

                        if not(round(cif_weight,8)==1.):
                                print("!!! cif weight=",round(cif_weight,8))
                                sys.exit()

                                for i in range(len(site_list)):
                                        print(site_list[i],valence_list[i])     #print site and valence
                                sys.exit()


                        if high_press_switch==True:
                                high_press_check+=cifname+"\n"

                        
                        analized_list+=cifname+"\n" #check analized cifname
                                

                        valencegraf=True
                        if valencegraf==True and valence_cut_point and valence_list:# create valence histgram
                                cation_valence_list=[]
                                vclass_width=0.1

                                for i in range(len(site_list)):
                                        if cation == site_list[i]:
                                                #print(site_list[i],valence_list[i])
                                                cation_valence_list.append(valence_list[i])

                                print(cation_valence_list)
                                
                                print(gm.distanse_divider(cation_valence_list[0],vclass_width))
                                print(len(vcontena))
                                #print(vcontena[gm.distanse_divider(cation_valence_list[0],vclass_width)])
                                for i in range(len(cation_valence_list)):
                                        if cation_valence_list[i]>10:
                                                mt.text_output("error unexpected valence"+str(cation_valence_list),cifname,"contena/errors")
                                                break
                                        vcontena[gm.distanse_divider(cation_valence_list[i],vclass_width)]  +=1/len(cation_valence_list)
                                        #except :
                                        #        print("!!!range error val=",cation_valence_list[i],gm.distanse_divider(cation_valence_list[i],class_width),max_distance/class_width)
                                        #        continue                
                                   

        """
        #print(high_press_check)
        #print(len(high_press_check.split("\n")))
        #print(cifname)
        mt.text_output(high_press_check,"high_press_Mg_1.40<distance<1.42","contena/errors")
        """

        mt.text_output(analized_list,'{f}_{M}_vcut={c}_bond={anion}_average={ave}_limsite{l}_Analized_list'.format(f=__file__,M=cation_list,c=valence_cut_point,O=nb_only_O,b=cn7_base,l=limit_site,anion=anion,ave=average_analize),"./")


        ### when cn>7 marge ce list
        more7marge=True
        if more7marge:
                print()
                print('marge')
                _contena={}
                for key in contena.keys():
                        _contena[key]={}
                        for ce in contena[key].keys():
                                cn=int(ce.split(":")[1])
                                if cn>=7:
                                        print(cn)
                                        if not str(cn) in _contena[key].keys(): 
                                                _contena[key][str(cn)]=contena[key][ce]
                                        else :
                                                for i in range(len(contena[key][ce])):
                                                        _contena[key][str(cn)][i]+=contena[key][ce][i]
                                else:
                                        _contena[key][ce]=contena[key][ce]  # cn=1~6
                #print(_contena['Ba2+'].keys())
                contena=dict(_contena)
                print()
        print("\nnew contena cn>7 ce data marge to cn label =\n",contena)

        ######################################################## Graph #####################################################################
        #axlist.append(plt.subplot(len(atomlist),1,atom_number+1)) 

        import seaborn as sns
        #sns.set_style("whitegrid")
        sns.set_style("darkgrid")

        base_colors=[(0.2,0.2,0.2),(0.75,0.75,0),(0.9,0,0.2),(0,0,1),(0.6,0.4,0.1),(0,0.4,0),
                     (0.3,0.6,1),(1,0.25,0.25),(0.4,0.5,1),(0.8,0.6,1),
                     (0,0.45,0.7),(0.6,1,0.1),(0,0,0),(0.5,0.5,0.5) ]
        #import shannon data
        with open('shannon1976radius') as l:
                shannon=l.read().split('\n')
        for i,line in enumerate(shannon[:len(shannon)]):
                shannon[i]=line.split()

        print("len contena",len(contena))

        axlist=[]           # list appended  for x asix
        #fig = plt.figure(figsize=(20,20/5*len(contena)))
        fig = plt.figure(figsize=(20+5.5,20/5*len(contena)))

        #fig.patch.set_facecolor('lightgray')    #backcolor

        newLabels, newHandles = [], []

        val_keys=[val for val in contena.keys()]
        print("val key=",val_keys)
        if valence_cut_point:
                #val_keys=gm.sort_ce(val_keys)
                val_keys=sorted(val_keys)
        print("val key=",val_keys)


        total_cif_count=0
        for v,cation in enumerate(val_keys):
                print("\n",cation)

                if fit_valence_list:
                        valence=cation
                else:
                        valence=""
                bottom=[0]*int(max_distance/class_width)       #for stack graf  **all stack

                # sort coordinaition numbar
                label=[]

                for cn in range(1,14):
                        for key in contena[cation].keys(): 
                                re_fd='.*:{d}\Z'.format(d=cn)
                                if re.findall(re_fd,key):
                                        label.append(key)
                                if key==str(cn) :
                                        label.append(key)
                print("label= ",label)
                label=gm.sort_ce(label)

                #color_set
                color=[gm.sym_to_color(c,cn7_base) for c in label]
                print(color)

                ### set graf ###
                #fig = plt.figure(figsize=(15,10))
                #fig.patch.set_facecolor('lightgray')    #backcolor
                plt.rcParams['axes.xmargin'] = 0        #(0,0)point
                ind=np.arange(int(max_distance/class_width))
                
                axlist.append(plt.subplot(len(contena.keys()),1,v+1))
                
                cif_count=0
                for i,lab in enumerate(label):
                                #plt.bar(ind,contena[cation][lab],color=color[i], label=label[i],bottom=bottom,width=1,align='edge' ,edgecolor=(0.25,0.25,0.25))   #create bar graf
                                #plt.bar(ind,contena[cation][lab],color=color[i], label=label[i],bottom=bottom,width=1,align='edge' ,edgecolor=(0.25,0.25,0.25),log=False,linewidth=0.8)   #create bar graf

                                #if ":" in lab:#nomal label ex "S:4"
                                if True:
                                        #plt.bar(ind,contena[cation][lab],color=color[i], label=label[i],bottom=bottom,width=1,align='edge' ,edgecolor=(0.25,0.25,0.25),log=False,linewidth=0.8)   #create bar graf
                                        plt.bar(ind,contena[cation][lab],color=color[i], label=label[i],bottom=bottom,width=1,align='edge' ,edgecolor=(0.,0.,0.),log=False,linewidth=1)   #create bar graf
                                        #plt.bar(ind,contena[cation][lab],color=color[i], label=label[i],bottom=bottom,width=1,align='edge' ,edgecolor=gm.inversion(color[i])  ,log=False,linewidth=1)   #create bar graf
                                else :        #marge file label ex "4"
                                        plt.bar(ind,contena[cation][lab],color=color[i], label=label[i],bottom=bottom,width=1,align='edge' ,edgecolor=(0,0,0),log=False,linewidth=0.8,hatch="//")   #create bar graf
                                        #plt.bar(ind,contena[cation][lab],color=color[i], label=label[i],bottom=bottom,width=1,align='edge' ,edgecolor=gm.inversion(color[i]) ,log=False,linewidth=0.8,hatch="xxx")   #create bar graf

                                bottom=np.array(bottom)+np.array(contena[cation][lab])                  #for stack graf
                                cif_count+=sum(contena[cation][lab])

                print("plot finish") 
                total_cif_count+=cif_count              
                    
                ### set graf labels                         
                #plt.legend(prop={'size':6,},title="coordination environment",loc='upper right', ncol=3,labelspacing=0,borderpad=0)

                if valence_cut_point:
                        #plt.title("{a} distance  {c}cif".format(a=cation,c=round(cif_count,2)),fontdict={'fontsize':24})
                        #plt.gca().set_title("{a} distance  {c}cif".format(a=cation,c=round(cif_count,2)), size=35,x=0.84,y=0.80,weight=400)
                        cation_str_sp=re.split(r"([A-Z][a-z]*)",cation)
                        print()
                        tex_formula=cation_str_sp[1]+"^{"+cation_str_sp[2]+"}"
                        cation_str="$\mathrm{"+tex_formula+"}$"
                        print(tex_formula)
                        print(cation_str)
                        plt.gca().set_title("{a}-{anion} length {c}cif".format(a=cation_str,c=round(cif_count),anion=anion), size=35,x=0.84,y=0.8,weight=400)
                else:
                        if not fit_valence:
                                title_valence=""
                        elif fit_valence[0]==1:
                                title_valence=""
                        else:
                                title_valence=fit_valence[0]



                        if not alcari:
                                #plt.title("{a} distance  {c}cif".format(a=cation,c=round(cif_count,2)),fontdict={'fontsize':24})
                                #plt.gca().set_title("{a} distance  {c}cif".format(a=cation,c=round(cif_count)), size=35,x=0.84,y=0.8,weight=400)
                                plt.gca().set_title("{a}-{anion} length {c}cif".format(a=cation,c=round(cif_count),anion=anion), size=35,x=0.84,y=0.8,weight=400)
                        else:
                                #plt.title("{a}{v}+ distance  {c}cif".format(a=cation,v=title_valence,c=int(cif_count)),fontdict={'fontsize':24})
                                #plt.gca().set_title("{a}{v}+ distance  {c}cif".format(a=cation,v=title_valence,c=int(cif_count)), size=35,x=0.84,y=0.8,weight=400)
                                #plt.gca().set_title("{a} distance  {c}cif".format(a=cation,c=round(cif_count)), size=35,x=0.84,y=0.8,weight=400)
                                plt.gca().set_title("{a}-{anion} length {c}cif".format(a=cation,c=round(cif_count),anion=anion), size=35,x=0.84,y=0.8,weight=400)

                if v==len(contena.keys())-1:            #most bottom label
                        #plt.xlabel("distance(Å)",fontdict={'fontsize':24},color=(0.3,0.4,0.5))
                        #plt.xlabel("Bond length (Å)",fontdict={'fontsize':24},color=(0.3,0.4,0.5))
                        plt.xlabel("Bond length (Å)",fontdict={'fontsize':30},color=(0.3,0.4,0.5))
                        axlist[v].xaxis.set_label_coords(0.5, -0.25)


                #plt.ylabel("weight",fontdict={'fontsize':24},color=(0.3,0.4,0.5))
                #plt.ylabel("Weight",fontdict={'fontsize':24},color=(0.3,0.4,0.5))
                plt.ylabel("Weight",fontdict={'fontsize':30},color=(0.3,0.4,0.5))
                #plt.tick_params(labelsize=18)
                label_point=np.arange(max_distance*10)/class_width/10
                #plt.xticks(label_point,np.arange(10*10)/10)

                if v==len(contena.keys())-1:            #most bottom label  
                        #plt.xticks(label_point,np.arange(10*10)/10)
                        plt.xticks(label_point,np.arange(10*10)/10,fontsize=30)
                else:
                        plt.xticks(label_point,[""]*(10*10))
                plt.tick_params(labelsize=30)
                plt.tight_layout()
                fig.align_labels()

                #print(np.log(np.arange(10*10)/10))
                
                ### plot shannon arrows
                """
                if valence=='':
                        print('valence=""')
                        for valence in fit_valence:
                                atom=cation +str(valence) +'+'
                                if not re.findall(r'[A-Z][a-z]*\d*\W',atom):continue
                                for line in shannon:
                                        if len(line)<=4:continue
                                        if re.findall(r'[A-Z][a-z]*\W',atom):   #for Li+
                                                atom=re.findall(r'[A-Z][a-z]*',atom)[0]+'1+'
                                        if (line[1]==re.findall(r'[A-Z][a-z]*',atom)[0]) and (re.findall(r'\d+',atom)[0] in line[2]):
                                                arrow_dict = dict(width=1,facecolor = base_colors[int(re.findall(r'\d+',line[3])[0])-1] ,edgecolor = base_colors[int(re.findall(r'\d+',line[3])[0])-1])
                                                text_dict = dict(boxstyle = "round",fc = 'white', ec = base_colors[int(re.findall(r'\d+',line[3])[0])-1],alpha=0.5)

                                                #axlist[v].annotate("{val}:val\n{cn}:cn".format(cn=line[3],val=line[2]), size = 10, color = "black",
                                                #                   xy = ((float(line[4])+R_O-0.14)/class_width, 0),
                                                #                   xytext = ((float(line[4])+R_O-0.14)/class_width+0.001,-max(bottom)/5),bbox = text_dict, arrowprops = arrow_dict)
                                                axlist[v].annotate("{val}:val\n{cn}:cn".format(cn=line[3],val=line[2]), size = 10, color = "black",
                                                                   xy = ((float(line[4])+R_O-0.14)/class_width, 0),
                                                                   xytext = ((float(line[4])+R_O-0.14)/class_width+0.001,-max(bottom)/4.5),bbox = text_dict, arrowprops = arrow_dict)
                                                
                                #print("max(bottom)",max(bottom))
                """
                ### ver.2
                ### valence label arrow
                arrow_dict = dict(width=1,facecolor = (0,0,0) ,edgecolor = (0,0,0),alpha=0)
                text_dict = dict(boxstyle = "round",fc = 'white', ec = (0,0,0),alpha=1)
                if not alcari and not fit_valence:
                        #arrow_text="      valence\ncoordination number"
                        pass
                else:
                        #arrow_text="coordination number"
                        pass

                #axlist[v].annotate(arrow_text, size = 13, color = "black",
                #                   xy = ((xlim[1]-(xlim[1]-xlim[0])/8.5)/class_width, 0),
                #                   xytext = ((xlim[1]-(xlim[1]-xlim[0])/8.51)/class_width,-max(bottom)/3.6),
                #                   bbox = text_dict, arrowprops = arrow_dict)
                #print(arrow_text)
                
                #if valence=='':                                 #when no valence ,plot all valence arrows.
                #arrow_textsize=13
                arrow_textsize=20
                if not valence_cut_point:
                        print('valence=""')

                        any_valence=np.arange(7)[1:]
                        arrow_list=[]
                        for valence in any_valence:
                                atom=cation +str(valence) +'+'
                                if not re.findall(r'[A-Z][a-z]*\d*\W',atom):continue
                                for line in shannon:
                                        if len(line)<=4:continue
                                        if re.findall(r'[A-Z][a-z]*\W',atom):   #for Li+
                                                atom=re.findall(r'[A-Z][a-z]*',atom)[0]+'1+'

                                        if (line[1]==re.findall(r'[A-Z][a-z]*',atom)[0]) and (re.findall(r'\d+',atom)[0] in line[2]):
                                                arrow_list.append(line)
                        print("arrow list",arrow_list)


                        ### plot arrows
                        for line in arrow_list:
                                arrow_dict = dict(width=1,facecolor = base_colors[int(re.findall(r'\d+',line[3])[0])-1] ,edgecolor = base_colors[int(re.findall(r'\d+',line[3])[0])-1])
                                text_dict = dict(boxstyle = "round",fc = 'white', ec = base_colors[int(re.findall(r'\d+',line[3])[0])-1],alpha=0.5)

                                #axlist[v].annotate("{val}:val\n{cn}:cn".format(cn=line[3],val=line[2]), size = 10, color = "black",
                                #                   xy = ((float(line[4])+R_O-0.14)/class_width, 0),
                                #                   xytext = ((float(line[4])+R_O-0.14)/class_width+0.001,-max(bottom)/5),bbox = text_dict, arrowprops = arrow_dict)

                                if not alcari:
                                        arrow_text="{val}:val\n{cn}:cn".format(cn=line[3],val=line[2])
                                else:         
                                        arrow_text="{cn}:cn".format(cn=line[3])

                                #if line==arrow_list[len(arrow_list)-1]:
                                #if line==arrow_list[0]:

                                
                                arrow_text=re.sub(r":[a-z]+",r"",arrow_text)    #arrow text  "3:val\n 6:cn" >> "3\n6"

                                if v==len(contena.keys())-1:            #most bottom label 
                                        axlist[v].annotate(arrow_text, size = arrow_textsize, color = "black",
                                                           xy = ((float(line[4])+R_O-0.14)/class_width, 0),
                                                           xytext = ((float(line[4])+R_O-0.14)/class_width+0.001,-max(bottom)/3.6),bbox = text_dict, arrowprops = arrow_dict),
                                                           #xytext = ((float(line[4])+R_O-0.14)/class_width+0.001,-max(bottom)/2),bbox = text_dict, arrowprops = arrow_dict)
                                else:
                                        axlist[v].annotate(arrow_text, size = arrow_textsize, color = "black",
                                                           xy = ((float(line[4])+R_O-0.14)/class_width, 0),
                                                           #xytext = ((float(line[4])+R_O-0.14)/class_width+0.001,-max(bottom)/6),bbox = text_dict, arrowprops = arrow_dict)
                                                           xytext = ((float(line[4])+R_O-0.14)/class_width+0.001,-max(bottom)/5),bbox = text_dict, arrowprops = arrow_dict)
                                                           

                                ##**** size=13 defolt
                                ##**** size=20 abst



                else:   ### valence True
                        atom=cation +str(valence) +'+'
                        if not re.findall(r'[A-Z][a-z]*\d*\W',atom):continue

                        ### plot arrows
                        for line in shannon:
                                if len(line)<=4:continue
                                if re.findall(r'[A-Z][a-z]*\W',atom):   #for Li+
                                        atom=re.findall(r'[A-Z][a-z]*',atom)[0]+'1+'
                                if (line[1]==re.findall(r'[A-Z][a-z]*',atom)[0]) and (re.findall(r'\d+',atom)[0] in line[2]):
                                        arrow_dict = dict(width=1,facecolor = base_colors[int(re.findall(r'\d+',line[3])[0])-1] ,edgecolor = base_colors[int(re.findall(r'\d+',line[3])[0])-1])
                                        text_dict = dict(boxstyle = "round",fc = 'white', ec = base_colors[int(re.findall(r'\d+',line[3])[0])-1],alpha=0.5)

                                        #axlist[v].annotate("{val}:val\n{cn}:cn".format(cn=line[3],val=line[2]), size = 10, color = "black",
                                        #                   xy = ((float(line[4])+R_O-0.14)/class_width, 0),
                                        #                   xytext = ((float(line[4])+R_O-0.14)/class_width+0.001,-max(bottom)/5),bbox = text_dict, arrowprops = arrow_dict)
                                        #axlist[v].annotate("cn\n{cn}".format(cn=line[3],val=line[2]), size = 10, color = "black",
                                        #                   xy = ((float(line[4])+R_O-0.14)/class_width, 0),
                                        #                   xytext = ((float(line[4])+R_O-0.14)/class_width+0.001,-max(bottom)/5),bbox = text_dict, arrowprops = arrow_dict)
                                        #arrow_text="{cn}:cn".format(cn=line[3])
                                        #arrow_text=re.sub(r":[a-z]+",r"",arrow_text)
                                        #axlist[v].annotate(arrow_text, size = 13, color = "black",
                                        #                   xy = ((float(line[4])+R_O-0.14)/class_width, 0),
                                        #                   xytext = ((float(line[4])+R_O-0.14)/class_width+0.001,-max(bottom)/3.6),bbox = text_dict, arrowprops = arrow_dict)

                                        if v==len(contena.keys())-1:            #most bottom label 
                                                arrow_text="{cn}:cn".format(cn=line[3])
                                                arrow_text=re.sub(r":[a-z]+",r"",arrow_text)
                                                axlist[v].annotate(arrow_text, size = arrow_textsize, color = "black",
                                                                   xy = ((float(line[4])+R_O-0.14)/class_width, 0),
                                                                   xytext = ((float(line[4])+R_O-0.14)/class_width+0.001,-max(bottom)/3.6),bbox = text_dict, arrowprops = arrow_dict)
                                                                   #xytext = ((float(line[4])+R_O-0.14)/class_width+0.001,-max(bottom)/2),bbox = text_dict, arrowprops = arrow_dict)
                                        
                                        else:            
                                                arrow_text="{cn}:cn".format(cn=line[3])
                                                arrow_text=re.sub(r":[a-z]+",r"",arrow_text)
                                                axlist[v].annotate(arrow_text, size = arrow_textsize, color = "black",
                                                                   xy = ((float(line[4])+R_O-0.14)/class_width, 0),
                                                                   #xytext = ((float(line[4])+R_O-0.14)/class_width+0.001,-max(bottom)/6),bbox = text_dict, arrowprops = arrow_dict)
                                                                   xytext = ((float(line[4])+R_O-0.14)/class_width+0.001,-max(bottom)/5),bbox = text_dict, arrowprops = arrow_dict)
                                                                   


                print('arrow plot finish')

                ### legend sorted ###      
                handles, labels = plt.gca().get_legend_handles_labels() #take a legend data
                for handle, label in zip(handles, labels):#sort legend
                  if label not in newLabels:
                    newLabels.append(label)
                    newHandles.append(handle)


        ### x asix pointing
        for ax in axlist:
                ax.set_xlim([xlim[0]/class_width,xlim[1]/class_width])

        #plt.subplots_adjust(hspace=0.5)
        #plt.subplots_adjust(hspace=1.5)
        #plt.subplots_adjust(hspace=0.4)
        plt.subplots_adjust(hspace=0.2)

        ### sort legend
        #print(newHandles)
        print("newLabels",newLabels)
        label_handles={}
        for i in range(len(newLabels)):
                label_handles[newLabels[i]]=newHandles[i]
        sorted_labels=gm.sort_ce(newLabels)
        sorted_handles=[label_handles[label] for label in sorted_labels]

        #reverse legend
        sorted_labels.reverse()
        sorted_handles.reverse()

        if not figset_name=="abst":
                #fig.legend(newHandles, newLabels,     prop={'size':15,},title="coordination environment",title_fontsize=15,loc='upper right',borderaxespad=8, ncol=3,labelspacing=0,borderpad=0)
                if alcari:
                        #fig.legend(sorted_handles,sorted_labels,     prop={'size':20,},title="coordination environment",title_fontsize=20,loc='lower left',borderaxespad=5, ncol=3,labelspacing=0,borderpad=0,framealpha=0.05,facecolor=(0.25,0.25,0.25))
                        #fig.legend(sorted_handles,sorted_labels,     prop={'size':25,},title="coordination environment",title_fontsize=20,loc='lower left',borderaxespad=5, ncol=1,labelspacing=0,borderpad=0,framealpha=1,facecolor=(0.9,0.9,0.9))
                        fig.legend(sorted_handles,sorted_labels,     prop={'size':25,},title="coordination\nenvironment",title_fontsize=30,loc='lower left',borderaxespad=5, ncol=1,labelspacing=0,borderpad=0,framealpha=1,facecolor=(0.9,0.9,0.9))
                else:
                        #fig.legend(sorted_handles,sorted_labels,     prop={'size':20,},title="coordination environment",title_fontsize=20,loc='upper left',borderaxespad=4, ncol=3,labelspacing=0,borderpad=0,framealpha=0.05,facecolor=(0.25,0.25,0.25))
                        #fig.legend(sorted_handles,sorted_labels,     prop={'size':25,},title="coordination environment",title_fontsize=20,loc='upper left',borderaxespad=4, ncol=1,labelspacing=0,borderpad=0,framealpha=1,facecolor=(0.9,0.9,0.9))
                        fig.legend(sorted_handles,sorted_labels,     prop={'size':25,},title="coordination\nenvironment",title_fontsize=30,loc='upper left',borderaxespad=4, ncol=1,labelspacing=0,borderpad=0,framealpha=1,facecolor=(0.9,0.9,0.9))
                #fig.suptitle("{c} distance".format(c=cation_list),fontdict={'fontsize':24})


        # addapt setting
        if figset_name=="abst":
                for iax in axlist:
                        iax.tick_params(labelsize=30)
                #fig.legend(sorted_handles,sorted_labels,     prop={'size':25,},title="coordination environment",title_fontsize=25,loc='upper left',borderaxespad=4, ncol=2,labelspacing=0,borderpad=0,framealpha=0.05,facecolor=(0.25,0.25,0.25))
                #fig.legend(sorted_handles,sorted_labels,     prop={'size':25,},title="coordination environment",title_fontsize=25,loc='upper left',borderaxespad=4, ncol=2,labelspacing=0,borderpad=0,framealpha=0.05,facecolor=(0.25,0.25,0.25))
                fig.legend(sorted_handles,sorted_labels,     prop={'size':30,},title="coordination\nenvironment",title_fontsize=25,loc='upper left',borderaxespad=4, ncol=2,labelspacing=0,borderpad=0,framealpha=0.05,facecolor=(0.25,0.25,0.25))

        #plt.tight_layout()

        ###save figure

        #create directory for save fig
        save_dirname="./figures"
        try:os.mkdir(save_dirname)
        except:pass

        save_dirname+=("/"+__file__)
        try:os.mkdir(save_dirname)
        except:pass 

        save_dirname+=("/"+anion)
        try:os.mkdir(save_dirname)
        except:pass       

        save_dirname+=("/"+"average={a}".format(a=average_analize))
        try:os.mkdir(save_dirname)
        except:pass
        
        #plt.savefig('hist_dist_ce_val_{M}.png'.format(M=cation),bbox_inches="tight", pad_inches=0.0 ,format="png")
        #plt.savefig('hist_dist_ce_2.6_{M}_nbO={O}_cn7color={b}.svg'.format(M=cation_list,c=valence_cut_point,O=nb_only_O,b=cn7_base), pad_inches=1.0 ,format="svg")
        #plt.savefig('hist_dist_ce_2.6_{M}_nbO={O}_cn7color={b}.png'.format(M=cation_list,c=valence_cut_point,O=nb_only_O,b=cn7_base), pad_inches=1.0 ,format="png")
        if value_log:
                plt.savefig('{f}_{M}_vcut={c}_log.png'.format(f=__file__,M=cation_list,c=valence_cut_point,O=nb_only_O,b=cn7_base), pad_inches=1.0 ,format="png")
                sys.exit()

        #plt.savefig('{f}_{M}_vcut={c}_bond={anion}_average={ave}_limsite{l}.svg'.format(f=__file__,M=cation_list,c=valence_cut_point,O=nb_only_O,b=cn7_base,l=limit_site,anion=anion,ave=average_analize), pad_inches=1.0 ,format="svg")
        #plt.savefig('{f}_{M}_vcut={c}_bond={anion}_average={ave}_limsite{l}.png'.format(f=__file__,M=cation_list,c=valence_cut_point,O=nb_only_O,b=cn7_base,l=limit_site,anion=anion,ave=average_analize), pad_inches=1.0 ,format="png")
        
        fig_name='{M}_vcut={c}_limsite{l}_{set}'.format(M=cation_list,c=valence_cut_point,l=limit_site,set=figset_name)
        plt.savefig(save_dirname+"/"+fig_name+".png", pad_inches=1.0,format="png")
        print("save image")

        ###remove dat file
        if True:
                import module_remover as mr
                mr.remove(r'\d+_out\.dat','./contena/VBA_dat')
                mr.remove(r'\d+_inp\.dat','./contena/VBA_dat')
                print("finish remove")

        # create valence graf
        if valencegraf and vcontena:
                import module_valencegraf as mvg
                mvg.vhist(vcontena,re.split(r'\d',cation)[0],int(total_cif_count),valence_cut_point)      
                print("save valence hist")  

        if log_scale:   # log scale
                for ax in axlist:
                        #ax.set_xlim([np.log(xlim[0]/class_width),np.log(xlim[1]/class_width)])
                        ax.set_xlim([1,xlim[1]/class_width])
                        #ax.set_xticks(np.linspace(0, np.pi * 4, 17), minor=True)
                        ax.grid(which="both",axis="x")
                        plt.xscale('log',basex=2)
                plt.savefig('hist_dist_ce_2.51_valR_bd={c}_{M}_nbO={O}_logscale.png'.format(M=cation,c=valence_cut_point,O=nb_only_O), pad_inches=0.0 ,format="png")

## main
anion_type = sys.argv[1]
anion_list=[anion_type]
#anion_list=["O","N","F"]
#figset_list=["alkali","alkali_arth","3d_1","3d_2","4d_1","4d_2"]
#figset_list=["4d_2"]
#figset_list=["abst"]
#figset_list=["valence"]
#figset_list=["alkali","alkali_arth","3d_1","3d_2","4d_1","4d_2","valence"]

if len(sys.argv) == 2:
	figset_list=["alkali","alkali_arth","3d_1","3d_2","4d_1","4d_2"]
else:
	figset_list=["valence"]

average_analize=False
for anion in anion_list:
        for figset_name in figset_list:
                make_bondhist(figset_name,anion,average_analize)
                
