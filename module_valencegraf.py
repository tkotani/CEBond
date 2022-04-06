import warnings
warnings.filterwarnings('ignore')   #except warnning!
import sys,os,subprocess
import collections
import numpy as np
import matplotlib.pyplot as plt
import re
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import module_valence_bond_analizer_RBVS as vba
import graf_module as gm

###  "keyword=aaa" take aaa
def keycut(content,keyword):   
         return content.split(keyword)[1].split('\n')[0]

### distance Discretization ###
def distanse_divider(distance):
        n=class_width
        i=0
        limit=class_width

        while round(limit,6)<=distance:
                #print(limit,distance)
                limit+=n
                i+=1
        return i

def vhist(contena,cation,cif_count,valence_cut_point=None):
        ### paramater for graf ###
        class_width=0.1
        max_distance=10
        xlim=[0,8]
        limit_site=500  #None or Number

        fig=plt.figure(figsize=(5*1.414,5))
        #fig.patch.set_facecolor('lightgray')    #backcolor
        import seaborn as sns
        sns.set_style("darkgrid")

        plt.rcParams['axes.xmargin'] = 0        #(0,0)point
        ind=np.arange(int(max_distance/class_width))
        ax=plt.subplot()

        plt.bar(ind,contena,width=1,align='edge' ,edgecolor=(0.25,0.25,0.25),color=(0.1,0.45,1))   #create bar graf
                 
        #graf label                         
        #plt.legend(prop={'size':6,})
        #plt.title("{ca} Valence  {c}cif \nfractional oqupancy site={f}".format(c=cif_count,ca=cation,f=fractional_occupancy_count),fontdict={'fontsize':24})
        plt.title("{ca} valence  {c}cif ".format(c=cif_count,ca=cation),fontdict={'fontsize':24})
        plt.xlabel("valence",fontdict={'fontsize':24})
        plt.ylabel("frequency",fontdict={'fontsize':24})
        plt.tick_params(labelsize=18)

        label_point=np.arange(max_distance*10)/class_width/10
        label_point=np.arange(max_distance*10)/class_width/1
        plt.xticks(label_point,np.arange(10*10)/1)

        plt.tight_layout()
        #arrow_dict = dict(width=1,facecolor =(1,0,0) ,edgecolor = (1,0,0))
        #text_dict  = dict(boxstyle = "round",fc = 'white', ec = (1,0,0))
        #ax.annotate("{r}Å".format(r=r), size = 14, color = "black",xy = (r/class_width, 0),xytext = (r/class_width, -max(hist_contena)/8),bbox = text_dict, arrowprops = arrow_dict)

        if valence_cut_point:
                for line_point in valence_cut_point:
                        plt.plot([line_point/class_width,line_point/class_width],[0, max(contena)*1.1], color=(0.8,0,0), linestyle='dashed')

        ax.set_xlim([xlim[0]/class_width,xlim[1]/class_width])
        #plt.subplots_adjust(hspace=0.5)

        #create directory for save fig
        save_dirname="./figures"
        try:os.mkdir(save_dirname)
        except:pass
        save_dirname+=("/valence_hist")
        try:os.mkdir(save_dirname)
        except:pass 
        
        fig_name='{M}_vcut={c}'.format(M=cation,c=valence_cut_point)
        plt.savefig(save_dirname+"/"+fig_name+".png", pad_inches=1,format="png")
        print("save image")
