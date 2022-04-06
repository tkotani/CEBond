#!/usr/bin/env python3
#Ver.2.0  remove high pressure Mg 2019-06-07
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
import logging 
from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import  SimplestChemenvStrategy, MultiWeightsChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments
import warnings
warnings.filterwarnings('ignore')   #except warnning!
import sys,os,subprocess


anion_type = sys.argv[1]
if anion_type == 'O':
        ### O paramater ###
        data="list_O.dat"
        anion=['O','O2-','O-','O0-','O2.033-','O0.333-','O0.667-','O0.5-']
        r=1.4
        #r=1.26

if anion_type == 'N':
        ### N paramater ###
        data='list_N.dat'
        anion=['N','N3-','N0+','N-','N+','N0.33-','N0.333-']
        r=1.46
        #r=1.32

if anion_type == 'F':
        ### F paramater ###
        data='list_F.dat' #_P-T'
        anion=['F','F1-']
        r=1.33
        #r=1.19

lattice_oqupancy=0.74 

#__________________________________________________________#

#import text saved chem data
with open(data) as l:
        contents=l.read()

contents=contents.split('\n')

### input skip list
import module_text_contena as mt
#import skiplist as sk
#skip_list= sk.skiplist

dirlist=[]
for line in contents[:len(contents)-1]:
        split=line.split()
        dirlist.append(split[0])

cwd = os.getcwd()
#basedir = '/home/sawada/cod/cif/' # ! cif base directory, CHENGE here!!
#basedir = '/home/iwamoto/lab_iwamoto/cod_analysis_iwamoto/COD/O/' #edit 20220318 iwamoto
basedir = './COD/'+sys.argv[1]
os.chdir(basedir)                 # move to basedirectry

#ss=open('setup_st_errors-test',mode='w')

# Setup the local geometry finder
lgf = LocalGeometryFinder()
lgf.setup_parameters(centering_type='centroid', include_central_site_in_centroid=True)
#Get the strategy from D. Waroquiers et al., Chem Mater., 2017, 29, 8346.
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import MultiWeightsChemenvStrategy
from pymatgen.io.cif import CifParser
#d
symbol_list=[]
import numpy as np
import matplotlib.pyplot as plt
import module_text_contena as mt

check_text=""

###  "keyword=aaa" take aaa
def keycut(content,keyword):   
         return content.split(keyword)[1].split('\n')[0]
#
#print(skip_list)
#sys.exit()
# for ciffile in dirlist:
#         cifname=ciffile.rsplit("/",1)[1].split(".")[0]
#         if cifname in skip_list:
#                 print("skip",cifname)
#                 continue
# sys.exit()

contena=[]
print('dirlist=',len(dirlist))
for ciffile in dirlist:
        #print('cif name=',ciffile)
        cifname=ciffile.rsplit("/",1)[1].split(".")[0]
        
        # if cifname in skip_list:
        #         print("skip",cifname)
        #         continue

        try:
                parser=CifParser(ciffile)
        except:
                #ff.write(ciffile)
                #sys.stdout.flush
                print(ciffile,' error_1')
                #continue
        try:
                struct=parser.get_structures()[0]
                #fff=struct.formula
                #aaa=struct.lattice.matrix
        except:
                #gg.write(ciffile)
                #sys.stdout.flush
                print(ciffile,' error_2')
                continue

        #print(struct.lattice.matrix)
        lattice_volume=np.linalg.det(struct.lattice.matrix)

        anion_count=0
        for isite in struct.sites:
                #print(isite.lattice)
                #print(isite._species._data.keys())
                #print(anion)
                #print(set(isite._species))
                for a in anion:
                        
                        if a in isite._species:
                                anion_count+=1

        if anion_count==0: 
                print("!!! anion None",ciffile,struct)
                print(isite._species)
                continue
        R=(lattice_volume*lattice_oqupancy/anion_count*3/4/np.pi)**(1/3)
        #print('R=',R)
        
        contena.append(R)
        print(cifname)
	
#        if R<1.05: 
#                #print('!!! radius ',ciffile,'\nR=',R,"\n",struct)
#                check_text=check_text+"\n\n"+ciffile+'\nR='+str(R)+"\n\n"+mt.text_input(ciffile)
      
os.chdir(cwd) #edit 20220318 iwamoto              
#mt.text_output(check_text,__file__+"_check") # R < 1.05  

fig = plt.figure(figsize=(20,5))
###paramater
class_width=0.01
max_distance=10
#xlim=[0,2.5]
xlim = [1.0,2.5] #edit 20210405 iwamoto

### distance Discretization
def distanse_divider(distance):
        n=class_width
        i=0
        while n<distance:
                n=n+class_width
                i=i+1
        return i

hist_contena=[0]*int(max_distance/class_width)
for d in contena:
        hist_contena[distanse_divider(d)]+=1
        
plt.rcParams['axes.xmargin'] = 0        #(0,0)point
ind=np.arange(int(max_distance/class_width))
ax=plt.subplot()
plt.bar(ind,hist_contena,width=1,align='edge' ,edgecolor=(0.25,0.25,0.25),color=(0.75,0.5,0.5))   #create bar graf
cif_count=len(contena)            
#graf label                         
#splt.legend(prop={'size':6,})
plt.title("{a} radius  {c}cif".format(c=cif_count,a=anion[0]),fontdict={'fontsize':24})
#plt.xlabel("distance(Å)",fontdict={'fontsize':24})
#plt.ylabel("latice count",fontdict={'fontsize':24})
plt.xlabel("R(Å)",fontdict={'fontsize':24})
plt.ylabel("Frequency",fontdict={'fontsize':24})
plt.tick_params(labelsize=18)
label_point=np.arange(max_distance*10)/class_width/10
plt.xticks(label_point,np.arange(10*10)/10)
#plt.tight_layout()

arrow_dict = dict(width=1,facecolor =(1,0,0) ,edgecolor = (1,0,0))
text_dict  = dict(boxstyle = "round",fc = 'white', ec = (1,0,0))
ax.annotate("{r}Å".format(r=r), size = 18, color = "black",xy = (r/class_width, 0),xytext = (r/class_width, -max(hist_contena)/8),bbox = text_dict, arrowprops = arrow_dict)

ax.set_xlim([xlim[0]/class_width,xlim[1]/class_width])
plt.subplots_adjust(hspace=0.5)

#os.chdir('/home/sawada/')
plt.savefig('{f}_{M}.svg'.format(M=anion[0],f=__file__),bbox_inches="tight", pad_inches=0.0 ,format="svg")
plt.savefig('{f}_{M}.png'.format(M=anion[0],f=__file__),bbox_inches="tight", pad_inches=0.0 ,format="png")
print("finish save fig")

