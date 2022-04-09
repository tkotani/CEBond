# Core part of BVS is from Prof.M.Kanzaki, Institute for Planetary Materials, Okayama University 
import math ,re ,sys,subprocess
import numpy as np
#import module_text_contena as mt
import os

###############################method1#########################################
##### ver.1
def BVS_Analizer_R(lattice_a,lattice_b,lattice_c,site_list,coord_list,R):
        bvs_list=[]
        global cation
        for i,isite in enumerate(site_list):
                #print("\nsite",i,isite)
                ivbs=0
                for j ,jsite in enumerate(site_list):
                        if i==j:
                                #print(isite,jsite)
                                continue
                        if isite==jsite:
                                #print(isite,jsite)
                                continue

                        if jsite=="O":
                                cation=isite
                                sign=1
                        elif isite=="O":
                                cation=jsite
                                sign=-1
                        else:
                                #print(isite,jsite)
                                continue

                        p12=[coord_list[j][0]-coord_list[i][0],
                             coord_list[j][1]-coord_list[i][1],
                             coord_list[j][2]-coord_list[i][2]]

                        #print("abs=",(p12[0]**2+p12[1]**2+p12[2]**2)**0.5)

                        vbs12=vb_analizer_distance_withinR(lattice_a,lattice_b,lattice_c,
                                                           p12,R,cation,sign)
                        ivbs+=vbs12

                bvs_list.append(ivbs)
        return bvs_list

def vb_analizer_distance_withinR(a,b,c,p12,R,cation,sign):
        # need "a.out" conpiled by gfortran from "neigbour1.F"
        
        # latice input
        #a,b,c,='0.5, 0.5, 1.0','0.5, 1.0, 0.5','1.0, 0.5, 0.5'

        #p12,R ='.5,.5,0'  # vectol isite to jsite 
        #R=     '1.0'      # R range

        #aaa=a +" " +b +" " +c +" " +p12 +" " +R
        aaa=str(a) +" " +str(b) +" " +str(c) +" " +str(p12) +" " +str(R)
        #
        #print(aaa)
        aaa=re.sub(r'\[|\]',r'',aaa)
        #print("aaa=",aaa)

        #os.system('echo '+aaa+' |./a.out' )
        distance_text=subprocess.getoutput('echo '+aaa+' |./a.out' )

        #print('a.out')
        #print(distance_text)
        #print(distance_text.split("\n")[:3])
        distance_list=distance_text.split("\n")[3:]
        #print(len(distance_list))

        n_dist_list=[float(d) for d in distance_list]
        #print("a.out=",n_dist_list)

        ro,b=take_vbparamater(cation)  # take  bv paramaters

        vbs12=0
        for r in n_dist_list:
                vb=sign *math.exp((ro-r)/b)
                vbs12+=vb
        return vbs12

##### Ver.2 ##############################################################
#    to fast
#import module_text_contena as mt 
def BVS_Analizer_R_2(cifname,R,site_list,lattice_a=None,lattice_b=None,lattice_c=None,coord_list=None):
        bvs_list=[]
        global cation

        path_w='{icif}_inp.dat'.format(icif=cifname)
        text=str(lattice_a)+str(lattice_b)+str(lattice_c)+"\n"+str(R)
        cation_dict={}
        sign_dict={}
        
        for i,isite in enumerate(site_list):
                #print("\nsite",i,isite)
                ivbs=0
                for j ,jsite in enumerate(site_list):
                        if i==j:#same site skip
                                #print(isite,jsite)
                                continue
                        if isite==jsite:#same atom skip
                                #print(isite,jsite)
                                continue
                        if jsite=="O":
                                cation=isite
                                sign=1
                        elif isite=="O":
                                cation=jsite
                                sign=-1
                        else:
                                #print(isite,jsite)
                                continue
                        if coord_list:
                                p12=[coord_list[j][0]-coord_list[i][0],
                                     coord_list[j][1]-coord_list[i][1],
                                     coord_list[j][2]-coord_list[i][2]]

                                text=text +"\n"+str(i)+" "+str(j)+" "+str(p12)

                        p12_lavel=str(i)+" "+str(j)
                        cation_dict[p12_lavel]=cation
                        sign_dict[p12_lavel]  =sign
        text=re.sub(r'\[|\]|\,',r' ',text)
        
        try:   
                with open("./VBA_dat/{icif}_out.dat".format(icif=cifname)) as f:
                        content=f.read()
                print("cash exist")

        except:
                save_dirname="./VBA_dat"
                try:os.mkdir(save_dirname)
                except:pass
                with open("./VBA_dat/"+path_w,mode='w') as f:
                        f.write(text)                          
                print("Fortran neib start")
                subprocess.getoutput('cd ./VBA_dat; echo {icif}|../a.out; cd ..'.format(icif=cifname))
                print("Fortran neib finish")
                with open("./VBA_dat/{icif}_out.dat".format(icif=cifname)) as f:
                        content=f.read()

        bvs_list=[0]*len(site_list)
        #print("content\n",content)
        for line in content.split('\n')[1:]:
                if line=='':continue
                sline=re.split(r'\s+',line)
                """
                for i in range(len(site_list)):
                        for j in range(len(site_list)):
                                if not sline[1]==str(i):continue
                                if not sline[2]==str(j):continue
                                for r in sline[3:]:
                                        r=float(r)
                                        cation=cation_dict[str(i)+" "+str(j)]
                                        sign  =sign_dict[str(i)+" "+str(j)]

                                        ro,b=take_vbparamater(cation) # take  bv paramaters
                                        vb=sign *math.exp((ro-r)/b)
                                        bvs_list[i]+=vb
                """
                
                cation=cation_dict[str(sline[1])+" "+str(sline[2])]
                sign  =sign_dict[str(sline[1])+" "+str(sline[2])]
                for r in sline[3:]:
                      r=float(r)
                      ro,b=take_vbparamater(cation) # take  bv paramaters
                      if ro==None or b==None:continue
                      vb=sign *math.exp((ro-r)/b)  
                      bvs_list[int(sline[1])]+=vb
        #print(bvs_list)
        #sys.exit()
        print("finish vb analize")
        return bvs_list

### Trance from cartegian to aaa
def trance_coordination_f_to_cartesian(a,b,c,alpha,beta,gamma):
        vec_a=[a,0,0] 
        vec_b=[b*np.cos(gamma/180*np.pi)  ,b*np.sin(gamma/180*np.pi)  ,0]
        vec_c=[c*np.cos(beta/180*np.pi),
               c*(np.cos(alpha/180.*np.pi)-np.cos(gamma/180*np.pi))  /np.sin(gamma/180*np.pi),
               c*(     np.sin(beta/180*np.pi)**2   -(   (np.cos(alpha/180.*np.pi)-np.cos(gamma/180*np.pi))  /np.sin(gamma/180*np.pi)  )**2    )**0.5]

        if (round((vec_a[0]**2+vec_a[1]**2+vec_a[2]**2)**0.5   -a,8)!=0)or(round((vec_b[0]**2+vec_b[1]**2+vec_b[2]**2)**0.5  -b,8)!=0) or(round((vec_c[0]**2+vec_c[1]**2+vec_c[2]**2)**0.5  -c,8)!=0):
                print("!!!Error !!! trance_coordination_f_to_cartesian")
                print(a,b,c)
                print((vec_a[0]**2+vec_a[1]**2+vec_a[2]**2)**0.5   ,(vec_b[0]**2+vec_b[1]**2+vec_b[2]**2)**0.5    ,(vec_c[0]**2+vec_c[1]**2+vec_c[2]**2)**0.5)
        return vec_a,vec_b,vec_c

#################################method2##########################################
## analize valence
## Search neibour in neighbor set
def vb_analize_from_site(center_symbol,center_position,neigbor_symbols,neighbor_positions):
        #print('<import vb_analize_from_site>')
        ### input ###
        #1.center site symbol(str object) and center position(list object)
        #2.neighbor site symbol(list object) and neighbor site position(list object)

        #print(center_symbol,center_position)

        center=re.match(r'[A-Z][a-z]*',center_symbol).group().strip() 
        neighbors=[re.match(r'[A-Z][a-z]*',symbol).group().strip() for symbol in neigbor_symbols] 
        #print(center,neighbors)

        ### get vb paramater for center site ###
        Ro_site,B_site  =  take_vbparamater(center_symbol)
        #print("Ro=",Ro_site)
        #print("B=",B_site)
        
        if (Ro_site==None) or (B_site==None):
                print("!! No vb paramater ") 
                return None

        ### get vb parameter for neighbor sites ###
        Ro_neib,B_neib=take_vbparamater_for_list(neigbor_symbols)
        #print("Ro_nb=",Ro_neib)
        #print("B_nb =",B_neib)

        if (None in set(Ro_neib)) or (None in set(B_neib)):
                print("!! No vb paramater ") 
                #sys.exit()
                return None

        vb_sum=0
        for i,neighbor in enumerate(neighbors):
                global ro,b
                if (center!=neighbor) and ((center=="O") or (neighbor=="O") ): #ex. cation-O or O-cation bond
                        
                        if center=="O":
                                ro=Ro_neib[i]
                                b =B_neib[i]
                                sign=-1
                        elif neighbor=='O':
                                ro=Ro_site
                                b=B_site
                                sign=1
                        else:
                                print("!!! unexpected bond")
                else:
                        print("!!!",center,"-",neighbor,'bond   bond convination skip vb analize')
                        continue
                print(center,"-",neighbor,"   ro,b=",ro,b)

                ### analize ###
                ### calcurate distance ###
                distance_x= neighbor_positions[i][0]-center_position[0]
                distance_y= neighbor_positions[i][1]-center_position[1]
                distance_z= neighbor_positions[i][2]-center_position[2]

                r=math.sqrt(distance_x**2 +distance_y**2 +distance_z**2)

                ### valence bond ###
                vb=sign *math.exp((ro-r)/b)
                print("vb=",vb)
                vb_sum+=vb

        round_vb_sum=round(vb_sum,1)
        return round_vb_sum


####################################paramater##################################
def take_vbparamater(site_symbol): 
        #print("<take_vb for site>")
        #print(site_symbol)
        #symbol=re.sub(r'\d',r'',site_symbol).strip()
        #kind=[symbol]
        
        kind=re.findall(r'[A-Z][a-z]*',site_symbol)

        R0=[None]
        B= [None]
        V= [None]

        for i in range(len(kind)):
                if kind[i] == 'H':
	                R0[i] = 0.918
	                B[i] = 0.427
	                V[i] = 1.0
                elif kind[i] == 'D': # same as H
	                R0[i] = 0.918
	                B[i] = 0.427
	                V[i] = 1.0
                elif kind[i] == 'Li':
	                R0[i] = 1.062
	                B[i] = 0.642
	                V[i] = 1.0
                elif kind[i] == 'Be':
	                R0[i] = 1.429
	                B[i] = 0.297
	                V[i] = 2.0
                elif kind[i] == 'B':
	                R0[i] = 1.372
	                B[i] = 0.357
	                V[i] = 3.0
                elif kind[i] == 'C':
	                R0[i] = 1.398
	                B[i] = 0.399
	                V[i] = 4.0
                elif kind[i] == 'N':
	                R0[i] = 1.492
	                B[i] = 0.482
	                V[i] = 5.0
                elif kind[i] == 'Na':
	                R0[i] = 1.695
	                B[i] = 0.420
	                V[i] = 1.0
                elif kind[i] == 'Mg':
	                R0[i] = 1.608
	                B[i] = 0.443
	                V[i] = 2.0
                elif kind[i] == 'Al':
	                R0[i] = 1.634
	                B[i] = 0.390
	                V[i] = 3.0
                elif kind[i] == 'Si':
	                R0[i] = 1.624
	                B[i] = 0.389
	                V[i] = 4.0
        #	elif kind[i] == 'P': # P3+
        #		polyvalent = 1
        #		R0[i] = 1.655
        #		B[i] = 0.399
        #		V[i] = 3.0
                elif kind[i] == 'P': # P5+
	                polyvalent = 1
	                R0[i] = 1.624
	                B[i] = 0.399
	                V[i] = 5.0
        #	elif kind[i] == 'S': # S4+
        #		polyvalent = 1
        #		R0[i] = 1.643
        #		B[i] = 0.399
        #		V[i] = 4.0
                elif kind[i] == 'S': # S6+
	                polyvalent = 1
	                R0[i] = 1.634
	                B[i] = 0.399
	                V[i] = 6.0
        #	elif kind[i] == 'Cl': # Cl3+
        #		polyvalent = 1
        #		R0[i] = 1.722
        #		B[i] = 0.370
        #		V[i] = 3.0
        #	elif kind[i] == 'Cl': # Cl5+
        #		polyvalent = 1
        #		R0[i] = 1.703
        #		B[i] = 0.428
        #		V[i] = 5.0
                elif kind[i] == 'Cl': # Cl7+
	                polyvalent = 1
	                R0[i] = 1.669
	                B[i] = 0.428
	                V[i] = 7.0
                elif kind[i] == 'K':
	                R0[i] = 2.047
	                B[i] = 0.398
	                V[i] = 1.0
                elif kind[i] == 'Ca':
	                R0[i] = 1.907
	                B[i] = 0.409
	                V[i] = 2.0
                elif kind[i] == 'Sc':
	                R0[i] = 1.780
	                B[i] = 0.452
	                V[i] = 3.0
        #	elif kind[i] == 'Ti': # Ti3+
        #		polyvalent = 1
        #		R0[i] = 1.654
        #		B[i] = 0.542
        #		V[i] = 3.0
                elif kind[i] == 'Ti': # Ti4+
	                polyvalent = 1
	                R0[i] = 1.819
	                B[i] = 0.342
	                V[i] = 4.0
                elif kind[i] == 'V': # V3+
	                polyvalent = 1
	                R0[i] = 1.718
	                B[i] = 0.412
	                V[i] = 3.0
        #	elif kind[i] == 'V': # V4+
        #		polyvalent = 1
        #		R0[i] = 1.776
        #		B[i] = 0.364
        #		V[i] = 4.0
        #	elif kind[i] == 'V': # V5+
        #		polyvalent = 1
        #		R0[i] = 1.799
        #		B[i] = 0.388
        #		V[i] = 5.0
        #	elif kind[i] == 'Cr': # Cr2+
        #		polyvalent = 1
        #		R0[i] = 1.761
        #		B[i] = 0.350
        #		V[i] = 2.0
                elif kind[i] == 'Cr': # Cr3+
	                polyvalent = 1
	                R0[i] = 1.725
	                B[i] = 0.361
	                V[i] = 3.0
        #	elif kind[i] == 'Cr': # Cr4+
        #		polyvalent = 1
        #		R0[i] = 1.783
        #		B[i] = 0.410
        #		V[i] = 4.0
        #	elif kind[i] == 'Cr': # Cr5+
        #		polyvalent = 1
        #		R0[i] = 1.777
        #		B[i] = 0.375
        #		V[i] = 5.0
        #	elif kind[i] == 'Cr': # Cr6+
        #		polyvalent = 1
        #		R0[i] = 1.799
        #		B[i] = 0.375
        #		V[i] = 6.0
                elif kind[i] == 'Mn': # Mn2+
	                polyvalent = 1
	                R0[i] = 1.740
	                B[i] = 0.417
	                V[i] = 2.0
        #	elif kind[i] == 'Mn': # Mn3+
        #		polyvalent = 1
        #		R0[i] = 1.823
        #		B[i] = 0.247
        #		V[i] = 3.0
        #	elif kind[i] == 'Mn': # Mn4+
        #		polyvalent = 1
        #		R0[i] = 1.750
        #		B[i] = 0.374
        #		V[i] = 4.0
        #	elif kind[i] == 'Mn': # Mn5+
        #		polyvalent = 1
        #		R0[i] = 1.781
        #		B[i] = 0.375
        #		V[i] = 5.0
        #	elif kind[i] == 'Mn': # Mn6+
        #		polyvalent = 1
        #		R0[i] = 1.814
        #		B[i] = 0.375
        #		V[i] = 6.0
        #	elif kind[i] == 'Mn': # Mn7+
        #		polyvalent = 1
        #		R0[i] = 1.819
        #		B[i] = 0.375
        #		V[i] = 7.0
                elif kind[i] == 'Fe': # Fe2+
	                polyvalent = 1
	                R0[i] = 1.658
	                B[i] = 0.447
	                V[i] = 2.0
        #	elif kind[i] == 'Fe': # Fe3+
        #		polyvalent = 1
        #		R0[i] = 1.766
        #		B[i] = 0.360
        #		V[i] = 3.0
                elif kind[i] == 'Co': # Co2+
	                polyvalent = 1
	                R0[i] = 1.698
	                B[i] = 0.376
	                V[i] = 2.0
        #	elif kind[i] == 'Co': # Co3+
        #		polyvalent = 1
        #		R0[i] = 1.655
        #		B[i] = 0.364
        #		V[i] = 3.0
        #	elif kind[i] == 'Co': # Co4+
        #		polyvalent = 1
        #		R0[i] = 1.729
        #		B[i] = 0.358
        #		V[i] = 4.0
                elif kind[i] == 'Ni': # Ni2+
	                polyvalent = 1
	                R0[i] = 1.689
	                B[i] = 0.347
	                V[i] = 2.0
        #	elif kind[i] == 'Ni': # Ni4+
        #		polyvalent = 1
        #		R0[i] = 1.734
        #		B[i] = 0.335
        #		V[i] = 4.0
        #	elif kind[i] == 'Cu': # Cu1+
        #		polyvalent = 1
        #		R0[i] = 1.601
        #		B[i] = 0.335
        #		V[i] = 1.0
                elif kind[i] == 'Cu': # Cu2+
	                polyvalent = 1
	                R0[i] = 1.687
	                B[i] = 0.355
	                V[i] = 2.0
        #	elif kind[i] == 'Cu': # Cu3+
        #		polyvalent = 1
        #		R0[i] = 1.737
        #		B[i] = 0.375
        #		V[i] = 3.0
                elif kind[i] == 'Zn':
	                R0[i] = 1.684
	                B[i] = 0.383
	                V[i] = 2.0
                elif kind[i] == 'Ga':
	                R0[i] = 1.736
	                B[i] = 0.345
	                V[i] = 3.0
                elif kind[i] == 'Ge':
	                R0[i] = 1.750
	                B[i] = 0.363
	                V[i] = 4.0
        #	elif kind[i] == 'As': # As3+
        #		polyvalent = 1
        #		R0[i] = 1.775
        #		B[i] = 0.423
        #		V[i] = 3.0
                elif kind[i] == 'As': # As5+
	                polyvalent = 1
	                R0[i] = 1.765
	                B[i] = 0.352
	                V[i] = 5.0
                elif kind[i] == 'Se': # Se4+
	                polyvalent = 1
	                R0[i] = 1.805
	                B[i] = 0.401
	                V[i] = 4.0
        #	elif kind[i] == 'Se': # Se6+
        #		polyvalent = 1
        #		R0[i] = 1.797
        #		B[i] = 0.399
        #		V[i] = 6.0
                elif kind[i] == 'Br': # Br5+ not anion
	                polyvalent = 1
	                R0[i] = 1.890
	                B[i] = 0.571
	                V[i] = 5.0
        #	elif kind[i] == 'Br': # Br7+ not anion
        #		polyvalent = 1
        #		R0[i] = 1.850
        #		B[i] = 0.428
        #		V[i] = 7.0
                elif kind[i] == 'Rb':
	                R0[i] = 1.993
	                B[i] = 0.478
	                V[i] = 1.0
                elif kind[i] == 'Sr':
	                R0[i] = 1.958
	                B[i] = 0.478
	                V[i] = 2.0
                elif kind[i] == 'Y':
	                R0[i] = 1.978
	                B[i] = 0.407
	                V[i] = 3.0
                elif kind[i] == 'Zr':
	                R0[i] = 1.913
	                B[i] = 0.406
	                V[i] = 4.0
        #	elif kind[i] == 'Nb': # Nb4+
        #		polyvalent = 1
        #		R0[i] = 1.853
        #		B[i] = 0.479
        #		V[i] = 4.0
                elif kind[i] == 'Nb': # Nb5+
	                polyvalent = 1
	                R0[i] = 1.909
	                B[i] = 0.369
	                V[i] = 5.0
        #	elif kind[i] == 'Mo': # Mo3+
        #		polyvalent = 1
        #		R0[i] = 1.792
        #		B[i] = 0.436
        #		V[i] = 3.0
        #	elif kind[i] == 'Mo': # Mo4+
        #		polyvalent = 1
        #		R0[i] = 1.834
        #		B[i] = 0.404
        #		V[i] = 4.0
        #	elif kind[i] == 'Mo': # Mo5+
        #		polyvalent = 1
        #		R0[i] = 1.888
        #		B[i] = 0.314
        #		V[i] = 5.0
                elif kind[i] == 'Mo': # Mo6+
	                polyvalent = 1
	                R0[i] = 1.903
	                B[i] = 0.349
	                V[i] = 6.0
                elif kind[i] == 'Tc':
	                R0[i] = 1.915
	                B[i] = 0.375
	                V[i] = 7.0
        #	elif kind[i] == 'Ru': # Ru3+
        #		polyvalent = 1
        #		R0[i] = 1.745
        #		B[i] = 0.401
        #		V[i] = 3.0
        #	elif kind[i] == 'Ru': # Ru4+
        #		polyvalent = 1
        #		R0[i] = 1.833
        #		B[i] = 0.366
        #		V[i] = 4.0
                elif kind[i] == 'Ru': # Ru5+
	                polyvalent = 1
	                R0[i] = 1.894
	                B[i] = 0.346
	                V[i] = 5.0
                elif kind[i] == 'Rh': # Rh3+
        #		polyvalent = 1
	                R0[i] = 1.769
	                B[i] = 0.369
	                V[i] = 3.0
        #	elif kind[i] == 'Rh': # Rh4+
        #		polyvalent = 1
        #		R0[i] = 1.836
        #		B[i] = 0.422
        #		V[i] = 4.0
                elif kind[i] == 'Pd': # Pd2+
	                polyvalent = 1
	                R0[i] = 1.749
	                B[i] = 0.375
	                V[i] = 2.0
        #	elif kind[i] == 'Pd': # Pd4+
        #		polyvalent = 1
        #		R0[i] = 1.856
        #		B[i] = 0.352
        #		V[i] = 4.0
                elif kind[i] == 'Ag':
	                R0[i] = 1.875
	                B[i] = 0.359
	                V[i] = 1.0
                elif kind[i] == 'Cd':
	                R0[i] = 1.827
	                B[i] = 0.430
	                V[i] = 2.0
                elif kind[i] == 'In':
	                R0[i] = 1.823
	                B[i] = 0.459
	                V[i] = 3.0
                elif kind[i] == 'Sn': # Sn2+
	                polyvalent = 1
	                R0[i] = 1.910
	                B[i] = 0.451
	                V[i] = 2.0
        #	elif kind[i] == 'Sn': # Sn4+
        #		polyvalent = 1
        #		R0[i] = 1.946
        #		B[i] = 0.274
        #		V[i] = 4.0
        #	elif kind[i] == 'Sb': # Sb3+
        #		polyvalent = 1
        #		R0[i] = 1.932
        #		B[i] = 0.435
        #		V[i] = 3.0
                elif kind[i] == 'Sb': # Sb5+
	                polyvalent = 1
	                R0[i] = 1.892
	                B[i] = 0.475
	                V[i] = 5.0
                elif kind[i] == 'Te': # Te4+
	                polyvalent = 1
	                R0[i] = 1.960
	                B[i] = 0.389
	                V[i] = 4.0
        #	elif kind[i] == 'Te': # Te6+
        #		polyvalent = 1
        #		R0[i] = 1.922
        #		B[i] = 0.387
        #		V[i] = 6.0
                elif kind[i] == 'I': # I5+
	                polyvalent = 1
	                R0[i] = 1.992
	                B[i] = 0.474
	                V[i] = 5.0
        #	elif kind[i] == 'I': # I7+
        #		polyvalent = 1
        #		R0[i] = 1.930
        #		B[i] = 0.299
        #		V[i] = 7.0
                elif kind[i] == 'Cs':
	                R0[i] = 2.296
	                B[i] = 0.411
	                V[i] = 1.0
                elif kind[i] == 'Ba':
	                R0[i] = 2.223
	                B[i] = 0.406
	                V[i] = 2.0
                elif kind[i] == 'La':
	                R0[i] = 2.179
	                B[i] = 0.359
	                V[i] = 3.0
                elif kind[i] == 'Ce': # Ce3+
	                polyvalent = 1
	                R0[i] = 2.114
	                B[i] = 0.389
	                V[i] = 3.0
        #	elif kind[i] == 'Ce': # Ce4+
        #		polyvalent = 1
        #		R0[i] = 2.046
        #		B[i] = 0.416
        #		V[i] = 4.0
                elif kind[i] == 'Pr':
	                R0[i] = 2.071
	                B[i] = 0.411
	                V[i] = 3.0
                elif kind[i] == 'Nd':
	                R0[i] = 2.103
	                B[i] = 0.371
	                V[i] = 3.0
                elif kind[i] == 'Sm':
	                R0[i] = 2.049
	                B[i] = 0.404
	                V[i] = 3.0
        #	elif kind[i] == 'Eu': # Eu2+
        #		polyvalent = 1
        #		R0[i] = 1.943
        #		B[i] = 0.490
        #		V[i] = 2.0
                elif kind[i] == 'Eu': # Eu3+
	                polyvalent = 1
	                R0[i] = 2.068
	                B[i] = 0.359
	                V[i] = 3.0
                elif kind[i] == 'Gd':
	                R0[i] = 1.988
	                B[i] = 0.433
	                V[i] = 3.0
                elif kind[i] == 'Tb': # Tb3+
	                polyvalent = 1
	                R0[i] = 2.020
	                B[i] = 0.379
	                V[i] = 3.0
        #	elif kind[i] == 'Tb': # Tb4+
        #		polyvalent = 1
        #		R0[i] = 2.018
        #		B[i] = 0.395
        #		V[i] = 4.0
                elif kind[i] == 'Dy':
	                R0[i] = 2.002
	                B[i] = 0.389
	                V[i] = 3.0
                elif kind[i] == 'Ho':
	                R0[i] = 1.993
	                B[i] = 0.387
	                V[i] = 3.0
                elif kind[i] == 'Er':
	                R0[i] = 1.991
	                B[i] = 0.373
	                V[i] = 3.0
                elif kind[i] == 'Tm':
	                R0[i] = 1.977
	                B[i] = 0.381
	                V[i] = 3.0
                elif kind[i] == 'Yb':
	                R0[i] = 1.969
	                B[i] = 0.373
	                V[i] = 3.0
                elif kind[i] == 'Lu':
	                R0[i] = 1.939
	                B[i] = 0.403
	                V[i] = 3.0
                elif kind[i] == 'Hf':
	                R0[i] = 1.923
	                B[i] = 0.375
	                V[i] = 4.0
                elif kind[i] == 'Ta':
	                R0[i] = 1.916
	                B[i] = 0.343
	                V[i] = 5.0
        #	elif kind[i] == 'W': # W5+
        #		polyvalent = 1
        #		R0[i] = 1.848
        #		B[i] = 0.553
        #		V[i] = 5.0
                elif kind[i] == 'W': # W6+
	                polyvalent = 1
	                R0[i] = 1.909
	                B[i] = 0.339
	                V[i] = 6.0
        #	elif kind[i] == 'Re': # Re5+
        #		polyvalent = 1
        #		R0[i] = 1.834
        #		B[i] = 0.557
        #		V[i] = 5.0
                elif kind[i] == 'Re': # Re7+
	                polyvalent = 1
	                R0[i] = 1.943
	                B[i] = 0.406
	                V[i] = 7.0
        #	elif kind[i] == 'Os': # Os5+
        #		polyvalent = 1
        #		R0[i] = 1.870
        #		B[i] = 0.485
        #		V[i] = 5.0
        #	elif kind[i] == 'Os': # Os6+
        #		polyvalent = 1
        #		R0[i] = 1.904
        #		B[i] = 0.375
        #		V[i] = 6.0
        #	elif kind[i] == 'Os': # Os7+
        #		polyvalent = 1
        #		R0[i] = 1.937
        #		B[i] = 0.349
        #		V[i] = 7.0
                elif kind[i] == 'Os': # Os8+
	                polyvalent = 1
	                R0[i] = 1.966
	                B[i] = 0.405
	                V[i] = 8.0
        #	elif kind[i] == 'Ir': # Ir3+
        #		polyvalent = 1
        #		R0[i] = 1.755
        #		B[i] = 0.414
        #		V[i] = 3.0
                elif kind[i] == 'Ir': # Ir4+
	                polyvalent = 1
	                R0[i] = 1.909
	                B[i] = 0.258
	                V[i] = 4.0
        #	elif kind[i] == 'Ir': # Ir5+
        #		polyvalent = 1
        #		R0[i] = 1.909
        #		B[i] = 0.449
        #		V[i] = 5.0
        #	elif kind[i] == 'Pt': # Pt2+
        #		polyvalent = 1
        #		R0[i] = 1.742
        #		B[i] = 0.375
        #		V[i] = 2.0
                elif kind[i] == 'Pt': # Pt4+
	                polyvalent = 1
	                R0[i] = 1.856
	                B[i] = 0.407
	                V[i] = 4.0
                elif kind[i] == 'Au': 
	                R0[i] = 1.890
	                B[i] = 0.375
	                V[i] = 3.0
                elif kind[i] == 'Hg': 
	                R0[i] = 1.947
	                B[i] = 0.370
	                V[i] = 2.0
                elif kind[i] == 'Tl': # Tl1+ 
	                polyvalent = 1
	                R0[i] = 2.063
	                B[i] = 0.422
	                V[i] = 1.0
        #	elif kind[i] == 'Tl': # Tl3+ 
        #		polyvalent = 1
        #		R0[i] = 1.874
        #		B[i] = 0.504
        #		V[i] = 3.0
                elif kind[i] == 'Pb': # Pb2+ 
	                polyvalent = 1
	                R0[i] = 2.032
	                B[i] = 0.442
	                V[i] = 2.0
        #	elif kind[i] == 'Pb': # Pb4+ 
        #		polyvalent = 1
        #		R0[i] = 2.056
        #		B[i] = 0.280
        #		V[i] = 4.0
                elif kind[i] == 'Bi': # Bi3+ 
	                polyvalent = 1
	                R0[i] = 2.068
	                B[i] = 0.389
	                V[i] = 3.0
        #	elif kind[i] == 'Bi': # Bi5+ 
        #		polyvalent = 1
        #		R0[i] = 2.050
        #		B[i] = 0.318
        #		V[i] = 5.0
                elif kind[i] == 'Th': 
	                R0[i] = 2.117
	                B[i] = 0.420
	                V[i] = 4.0
        #	elif kind[i] == 'U': # U4+ 
        #		polyvalent = 1
        #		R0[i] = 2.100
        #		B[i] = 0.373
        #		V[i] = 4.0
        #	elif kind[i] == 'U': # U5+ 
        #		polyvalent = 1
        #		R0[i] = 2.009
        #		B[i] = 0.660
        #		V[i] = 5.0
                elif kind[i] == 'U': # U6+ 
	                polyvalent = 1
	                R0[i] = 2.046
	                B[i] = 0.473
	                V[i] = 6.0
                elif kind[i] == 'Np': # Np5+ 
	                polyvalent = 1
	                R0[i] = 2.036
	                B[i] = 0.411
	                V[i] = 5.0
        #	elif kind[i] == 'Np': # Np6+ 
        #		polyvalent = 1
        #		R0[i] = 2.022
        #		B[i] = 0.523
        #		V[i] = 6.0
        #	elif kind[i] == 'Np': # Np7+ 
        #		polyvalent = 1
        #		R0[i] = 2.076
        #		B[i] = 0.477
        #		V[i] = 7.0
                elif kind[i] == 'Am': 
	                R0[i] = 2.068
	                B[i] = 0.392
	                V[i] = 3.0
                elif kind[i] == 'Cm': 
	                R0[i] = 2.034
	                B[i] = 0.412
	                V[i] = 3.0
                elif kind[i] == 'O':  # as O is anion
	                R0[i] = 0.0
	                B[i] = 0.0
	                V[i] = 2.0
        # not found in list
                else:
                        print('Not found in build-in list!  kind= ',kind[i])

                #if kind[i] != 'O':
                        #print(kind[i] + ': R0 = ' + str(R0[i]) + ' : B = ' + str(B[i]) + ' : Q = ' + str(V[i]))


                
        return R0[0] ,B[0] 


def take_vbparamater_for_list(site_list):
        #print("<take_vb_for nv>")
        #print(site_list)

        kind=[re.match(r'[A-Z][a-z]*',symbol).group().strip() for symbol in site_list]
        #print(kind)
        #sys.exit()
        R0=[None for i in range(len(site_list))]
        B= [None for i in range(len(site_list))]
        V= [None for i in range(len(site_list))]


        for i in range(len(kind)):
	        if kind[i] == 'H':
		        R0[i] = 0.918
		        B[i] = 0.427
		        V[i] = 1.0
	        elif kind[i] == 'D': # same as H
		        R0[i] = 0.918
		        B[i] = 0.427
		        V[i] = 1.0
	        elif kind[i] == 'Li':
		        R0[i] = 1.062
		        B[i] = 0.642
		        V[i] = 1.0
	        elif kind[i] == 'Be':
		        R0[i] = 1.429
		        B[i] = 0.297
		        V[i] = 2.0
	        elif kind[i] == 'B':
		        R0[i] = 1.372
		        B[i] = 0.357
		        V[i] = 3.0
	        elif kind[i] == 'C':
		        R0[i] = 1.398
		        B[i] = 0.399
		        V[i] = 4.0
	        elif kind[i] == 'N':
		        R0[i] = 1.492
		        B[i] = 0.482
		        V[i] = 5.0
	        elif kind[i] == 'Na':
		        R0[i] = 1.695
		        B[i] = 0.420
		        V[i] = 1.0
	        elif kind[i] == 'Mg':
		        R0[i] = 1.608
		        B[i] = 0.443
		        V[i] = 2.0
	        elif kind[i] == 'Al':
		        R0[i] = 1.634
		        B[i] = 0.390
		        V[i] = 3.0
	        elif kind[i] == 'Si':
		        R0[i] = 1.624
		        B[i] = 0.389
		        V[i] = 4.0
        #	elif kind[i] == 'P': # P3+
        #		polyvalent = 1
        #		R0[i] = 1.655
        #		B[i] = 0.399
        #		V[i] = 3.0
	        elif kind[i] == 'P': # P5+
		        polyvalent = 1
		        R0[i] = 1.624
		        B[i] = 0.399
		        V[i] = 5.0
        #	elif kind[i] == 'S': # S4+
        #		polyvalent = 1
        #		R0[i] = 1.643
        #		B[i] = 0.399
        #		V[i] = 4.0
	        elif kind[i] == 'S': # S6+
		        polyvalent = 1
		        R0[i] = 1.634
		        B[i] = 0.399
		        V[i] = 6.0
        #	elif kind[i] == 'Cl': # Cl3+
        #		polyvalent = 1
        #		R0[i] = 1.722
        #		B[i] = 0.370
        #		V[i] = 3.0
        #	elif kind[i] == 'Cl': # Cl5+
        #		polyvalent = 1
        #		R0[i] = 1.703
        #		B[i] = 0.428
        #		V[i] = 5.0
	        elif kind[i] == 'Cl': # Cl7+
		        polyvalent = 1
		        R0[i] = 1.669
		        B[i] = 0.428
		        V[i] = 7.0
	        elif kind[i] == 'K':
		        R0[i] = 2.047
		        B[i] = 0.398
		        V[i] = 1.0
	        elif kind[i] == 'Ca':
		        R0[i] = 1.907
		        B[i] = 0.409
		        V[i] = 2.0
	        elif kind[i] == 'Sc':
		        R0[i] = 1.780
		        B[i] = 0.452
		        V[i] = 3.0
        #	elif kind[i] == 'Ti': # Ti3+
        #		polyvalent = 1
        #		R0[i] = 1.654
        #		B[i] = 0.542
        #		V[i] = 3.0
	        elif kind[i] == 'Ti': # Ti4+
		        polyvalent = 1
		        R0[i] = 1.819
		        B[i] = 0.342
		        V[i] = 4.0
	        elif kind[i] == 'V': # V3+
		        polyvalent = 1
		        R0[i] = 1.718
		        B[i] = 0.412
		        V[i] = 3.0
        #	elif kind[i] == 'V': # V4+
        #		polyvalent = 1
        #		R0[i] = 1.776
        #		B[i] = 0.364
        #		V[i] = 4.0
        #	elif kind[i] == 'V': # V5+
        #		polyvalent = 1
        #		R0[i] = 1.799
        #		B[i] = 0.388
        #		V[i] = 5.0
        #	elif kind[i] == 'Cr': # Cr2+
        #		polyvalent = 1
        #		R0[i] = 1.761
        #		B[i] = 0.350
        #		V[i] = 2.0
	        elif kind[i] == 'Cr': # Cr3+
		        polyvalent = 1
		        R0[i] = 1.725
		        B[i] = 0.361
		        V[i] = 3.0
        #	elif kind[i] == 'Cr': # Cr4+
        #		polyvalent = 1
        #		R0[i] = 1.783
        #		B[i] = 0.410
        #		V[i] = 4.0
        #	elif kind[i] == 'Cr': # Cr5+
        #		polyvalent = 1
        #		R0[i] = 1.777
        #		B[i] = 0.375
        #		V[i] = 5.0
        #	elif kind[i] == 'Cr': # Cr6+
        #		polyvalent = 1
        #		R0[i] = 1.799
        #		B[i] = 0.375
        #		V[i] = 6.0
	        elif kind[i] == 'Mn': # Mn2+
		        polyvalent = 1
		        R0[i] = 1.740
		        B[i] = 0.417
		        V[i] = 2.0
        #	elif kind[i] == 'Mn': # Mn3+
        #		polyvalent = 1
        #		R0[i] = 1.823
        #		B[i] = 0.247
        #		V[i] = 3.0
        #	elif kind[i] == 'Mn': # Mn4+
        #		polyvalent = 1
        #		R0[i] = 1.750
        #		B[i] = 0.374
        #		V[i] = 4.0
        #	elif kind[i] == 'Mn': # Mn5+
        #		polyvalent = 1
        #		R0[i] = 1.781
        #		B[i] = 0.375
        #		V[i] = 5.0
        #	elif kind[i] == 'Mn': # Mn6+
        #		polyvalent = 1
        #		R0[i] = 1.814
        #		B[i] = 0.375
        #		V[i] = 6.0
        #	elif kind[i] == 'Mn': # Mn7+
        #		polyvalent = 1
        #		R0[i] = 1.819
        #		B[i] = 0.375
        #		V[i] = 7.0
	        elif kind[i] == 'Fe': # Fe2+
		        polyvalent = 1
		        R0[i] = 1.658
		        B[i] = 0.447
		        V[i] = 2.0
        #	elif kind[i] == 'Fe': # Fe3+
        #		polyvalent = 1
        #		R0[i] = 1.766
        #		B[i] = 0.360
        #		V[i] = 3.0
	        elif kind[i] == 'Co': # Co2+
		        polyvalent = 1
		        R0[i] = 1.698
		        B[i] = 0.376
		        V[i] = 2.0
        #	elif kind[i] == 'Co': # Co3+
        #		polyvalent = 1
        #		R0[i] = 1.655
        #		B[i] = 0.364
        #		V[i] = 3.0
        #	elif kind[i] == 'Co': # Co4+
        #		polyvalent = 1
        #		R0[i] = 1.729
        #		B[i] = 0.358
        #		V[i] = 4.0
	        elif kind[i] == 'Ni': # Ni2+
		        polyvalent = 1
		        R0[i] = 1.689
		        B[i] = 0.347
		        V[i] = 2.0
        #	elif kind[i] == 'Ni': # Ni4+
        #		polyvalent = 1
        #		R0[i] = 1.734
        #		B[i] = 0.335
        #		V[i] = 4.0
        #	elif kind[i] == 'Cu': # Cu1+
        #		polyvalent = 1
        #		R0[i] = 1.601
        #		B[i] = 0.335
        #		V[i] = 1.0
	        elif kind[i] == 'Cu': # Cu2+
		        polyvalent = 1
		        R0[i] = 1.687
		        B[i] = 0.355
		        V[i] = 2.0
        #	elif kind[i] == 'Cu': # Cu3+
        #		polyvalent = 1
        #		R0[i] = 1.737
        #		B[i] = 0.375
        #		V[i] = 3.0
	        elif kind[i] == 'Zn':
		        R0[i] = 1.684
		        B[i] = 0.383
		        V[i] = 2.0
	        elif kind[i] == 'Ga':
		        R0[i] = 1.736
		        B[i] = 0.345
		        V[i] = 3.0
	        elif kind[i] == 'Ge':
		        R0[i] = 1.750
		        B[i] = 0.363
		        V[i] = 4.0
        #	elif kind[i] == 'As': # As3+
        #		polyvalent = 1
        #		R0[i] = 1.775
        #		B[i] = 0.423
        #		V[i] = 3.0
	        elif kind[i] == 'As': # As5+
		        polyvalent = 1
		        R0[i] = 1.765
		        B[i] = 0.352
		        V[i] = 5.0
	        elif kind[i] == 'Se': # Se4+
		        polyvalent = 1
		        R0[i] = 1.805
		        B[i] = 0.401
		        V[i] = 4.0
        #	elif kind[i] == 'Se': # Se6+
        #		polyvalent = 1
        #		R0[i] = 1.797
        #		B[i] = 0.399
        #		V[i] = 6.0
	        elif kind[i] == 'Br': # Br5+ not anion
		        polyvalent = 1
		        R0[i] = 1.890
		        B[i] = 0.571
		        V[i] = 5.0
        #	elif kind[i] == 'Br': # Br7+ not anion
        #		polyvalent = 1
        #		R0[i] = 1.850
        #		B[i] = 0.428
        #		V[i] = 7.0
	        elif kind[i] == 'Rb':
		        R0[i] = 1.993
		        B[i] = 0.478
		        V[i] = 1.0
	        elif kind[i] == 'Sr':
		        R0[i] = 1.958
		        B[i] = 0.478
		        V[i] = 2.0
	        elif kind[i] == 'Y':
		        R0[i] = 1.978
		        B[i] = 0.407
		        V[i] = 3.0
	        elif kind[i] == 'Zr':
		        R0[i] = 1.913
		        B[i] = 0.406
		        V[i] = 4.0
        #	elif kind[i] == 'Nb': # Nb4+
        #		polyvalent = 1
        #		R0[i] = 1.853
        #		B[i] = 0.479
        #		V[i] = 4.0
	        elif kind[i] == 'Nb': # Nb5+
		        polyvalent = 1
		        R0[i] = 1.909
		        B[i] = 0.369
		        V[i] = 5.0
        #	elif kind[i] == 'Mo': # Mo3+
        #		polyvalent = 1
        #		R0[i] = 1.792
        #		B[i] = 0.436
        #		V[i] = 3.0
        #	elif kind[i] == 'Mo': # Mo4+
        #		polyvalent = 1
        #		R0[i] = 1.834
        #		B[i] = 0.404
        #		V[i] = 4.0
        #	elif kind[i] == 'Mo': # Mo5+
        #		polyvalent = 1
        #		R0[i] = 1.888
        #		B[i] = 0.314
        #		V[i] = 5.0
	        elif kind[i] == 'Mo': # Mo6+
		        polyvalent = 1
		        R0[i] = 1.903
		        B[i] = 0.349
		        V[i] = 6.0
	        elif kind[i] == 'Tc':
		        R0[i] = 1.915
		        B[i] = 0.375
		        V[i] = 7.0
        #	elif kind[i] == 'Ru': # Ru3+
        #		polyvalent = 1
        #		R0[i] = 1.745
        #		B[i] = 0.401
        #		V[i] = 3.0
        #	elif kind[i] == 'Ru': # Ru4+
        #		polyvalent = 1
        #		R0[i] = 1.833
        #		B[i] = 0.366
        #		V[i] = 4.0
	        elif kind[i] == 'Ru': # Ru5+
		        polyvalent = 1
		        R0[i] = 1.894
		        B[i] = 0.346
		        V[i] = 5.0
	        elif kind[i] == 'Rh': # Rh3+
        #		polyvalent = 1
		        R0[i] = 1.769
		        B[i] = 0.369
		        V[i] = 3.0
        #	elif kind[i] == 'Rh': # Rh4+
        #		polyvalent = 1
        #		R0[i] = 1.836
        #		B[i] = 0.422
        #		V[i] = 4.0
	        elif kind[i] == 'Pd': # Pd2+
		        polyvalent = 1
		        R0[i] = 1.749
		        B[i] = 0.375
		        V[i] = 2.0
        #	elif kind[i] == 'Pd': # Pd4+
        #		polyvalent = 1
        #		R0[i] = 1.856
        #		B[i] = 0.352
        #		V[i] = 4.0
	        elif kind[i] == 'Ag':
		        R0[i] = 1.875
		        B[i] = 0.359
		        V[i] = 1.0
	        elif kind[i] == 'Cd':
		        R0[i] = 1.827
		        B[i] = 0.430
		        V[i] = 2.0
	        elif kind[i] == 'In':
		        R0[i] = 1.823
		        B[i] = 0.459
		        V[i] = 3.0
	        elif kind[i] == 'Sn': # Sn2+
		        polyvalent = 1
		        R0[i] = 1.910
		        B[i] = 0.451
		        V[i] = 2.0
        #	elif kind[i] == 'Sn': # Sn4+
        #		polyvalent = 1
        #		R0[i] = 1.946
        #		B[i] = 0.274
        #		V[i] = 4.0
        #	elif kind[i] == 'Sb': # Sb3+
        #		polyvalent = 1
        #		R0[i] = 1.932
        #		B[i] = 0.435
        #		V[i] = 3.0
	        elif kind[i] == 'Sb': # Sb5+
		        polyvalent = 1
		        R0[i] = 1.892
		        B[i] = 0.475
		        V[i] = 5.0
	        elif kind[i] == 'Te': # Te4+
		        polyvalent = 1
		        R0[i] = 1.960
		        B[i] = 0.389
		        V[i] = 4.0
        #	elif kind[i] == 'Te': # Te6+
        #		polyvalent = 1
        #		R0[i] = 1.922
        #		B[i] = 0.387
        #		V[i] = 6.0
	        elif kind[i] == 'I': # I5+
		        polyvalent = 1
		        R0[i] = 1.992
		        B[i] = 0.474
		        V[i] = 5.0
        #	elif kind[i] == 'I': # I7+
        #		polyvalent = 1
        #		R0[i] = 1.930
        #		B[i] = 0.299
        #		V[i] = 7.0
	        elif kind[i] == 'Cs':
		        R0[i] = 2.296
		        B[i] = 0.411
		        V[i] = 1.0
	        elif kind[i] == 'Ba':
		        R0[i] = 2.223
		        B[i] = 0.406
		        V[i] = 2.0
	        elif kind[i] == 'La':
		        R0[i] = 2.179
		        B[i] = 0.359
		        V[i] = 3.0
	        elif kind[i] == 'Ce': # Ce3+
		        polyvalent = 1
		        R0[i] = 2.114
		        B[i] = 0.389
		        V[i] = 3.0
        #	elif kind[i] == 'Ce': # Ce4+
        #		polyvalent = 1
        #		R0[i] = 2.046
        #		B[i] = 0.416
        #		V[i] = 4.0
	        elif kind[i] == 'Pr':
		        R0[i] = 2.071
		        B[i] = 0.411
		        V[i] = 3.0
	        elif kind[i] == 'Nd':
		        R0[i] = 2.103
		        B[i] = 0.371
		        V[i] = 3.0
	        elif kind[i] == 'Sm':
		        R0[i] = 2.049
		        B[i] = 0.404
		        V[i] = 3.0
        #	elif kind[i] == 'Eu': # Eu2+
        #		polyvalent = 1
        #		R0[i] = 1.943
        #		B[i] = 0.490
        #		V[i] = 2.0
	        elif kind[i] == 'Eu': # Eu3+
		        polyvalent = 1
		        R0[i] = 2.068
		        B[i] = 0.359
		        V[i] = 3.0
	        elif kind[i] == 'Gd':
		        R0[i] = 1.988
		        B[i] = 0.433
		        V[i] = 3.0
	        elif kind[i] == 'Tb': # Tb3+
		        polyvalent = 1
		        R0[i] = 2.020
		        B[i] = 0.379
		        V[i] = 3.0
        #	elif kind[i] == 'Tb': # Tb4+
        #		polyvalent = 1
        #		R0[i] = 2.018
        #		B[i] = 0.395
        #		V[i] = 4.0
	        elif kind[i] == 'Dy':
		        R0[i] = 2.002
		        B[i] = 0.389
		        V[i] = 3.0
	        elif kind[i] == 'Ho':
		        R0[i] = 1.993
		        B[i] = 0.387
		        V[i] = 3.0
	        elif kind[i] == 'Er':
		        R0[i] = 1.991
		        B[i] = 0.373
		        V[i] = 3.0
	        elif kind[i] == 'Tm':
		        R0[i] = 1.977
		        B[i] = 0.381
		        V[i] = 3.0
	        elif kind[i] == 'Yb':
		        R0[i] = 1.969
		        B[i] = 0.373
		        V[i] = 3.0
	        elif kind[i] == 'Lu':
		        R0[i] = 1.939
		        B[i] = 0.403
		        V[i] = 3.0
	        elif kind[i] == 'Hf':
		        R0[i] = 1.923
		        B[i] = 0.375
		        V[i] = 4.0
	        elif kind[i] == 'Ta':
		        R0[i] = 1.916
		        B[i] = 0.343
		        V[i] = 5.0
        #	elif kind[i] == 'W': # W5+
        #		polyvalent = 1
        #		R0[i] = 1.848
        #		B[i] = 0.553
        #		V[i] = 5.0
	        elif kind[i] == 'W': # W6+
		        polyvalent = 1
		        R0[i] = 1.909
		        B[i] = 0.339
		        V[i] = 6.0
        #	elif kind[i] == 'Re': # Re5+
        #		polyvalent = 1
        #		R0[i] = 1.834
        #		B[i] = 0.557
        #		V[i] = 5.0
	        elif kind[i] == 'Re': # Re7+
		        polyvalent = 1
		        R0[i] = 1.943
		        B[i] = 0.406
		        V[i] = 7.0
        #	elif kind[i] == 'Os': # Os5+
        #		polyvalent = 1
        #		R0[i] = 1.870
        #		B[i] = 0.485
        #		V[i] = 5.0
        #	elif kind[i] == 'Os': # Os6+
        #		polyvalent = 1
        #		R0[i] = 1.904
        #		B[i] = 0.375
        #		V[i] = 6.0
        #	elif kind[i] == 'Os': # Os7+
        #		polyvalent = 1
        #		R0[i] = 1.937
        #		B[i] = 0.349
        #		V[i] = 7.0
	        elif kind[i] == 'Os': # Os8+
		        polyvalent = 1
		        R0[i] = 1.966
		        B[i] = 0.405
		        V[i] = 8.0
        #	elif kind[i] == 'Ir': # Ir3+
        #		polyvalent = 1
        #		R0[i] = 1.755
        #		B[i] = 0.414
        #		V[i] = 3.0
	        elif kind[i] == 'Ir': # Ir4+
		        polyvalent = 1
		        R0[i] = 1.909
		        B[i] = 0.258
		        V[i] = 4.0
        #	elif kind[i] == 'Ir': # Ir5+
        #		polyvalent = 1
        #		R0[i] = 1.909
        #		B[i] = 0.449
        #		V[i] = 5.0
        #	elif kind[i] == 'Pt': # Pt2+
        #		polyvalent = 1
        #		R0[i] = 1.742
        #		B[i] = 0.375
        #		V[i] = 2.0
	        elif kind[i] == 'Pt': # Pt4+
		        polyvalent = 1
		        R0[i] = 1.856
		        B[i] = 0.407
		        V[i] = 4.0
	        elif kind[i] == 'Au': 
		        R0[i] = 1.890
		        B[i] = 0.375
		        V[i] = 3.0
	        elif kind[i] == 'Hg': 
		        R0[i] = 1.947
		        B[i] = 0.370
		        V[i] = 2.0
	        elif kind[i] == 'Tl': # Tl1+ 
		        polyvalent = 1
		        R0[i] = 2.063
		        B[i] = 0.422
		        V[i] = 1.0
        #	elif kind[i] == 'Tl': # Tl3+ 
        #		polyvalent = 1
        #		R0[i] = 1.874
        #		B[i] = 0.504
        #		V[i] = 3.0
	        elif kind[i] == 'Pb': # Pb2+ 
		        polyvalent = 1
		        R0[i] = 2.032
		        B[i] = 0.442
		        V[i] = 2.0
        #	elif kind[i] == 'Pb': # Pb4+ 
        #		polyvalent = 1
        #		R0[i] = 2.056
        #		B[i] = 0.280
        #		V[i] = 4.0
	        elif kind[i] == 'Bi': # Bi3+ 
		        polyvalent = 1
		        R0[i] = 2.068
		        B[i] = 0.389
		        V[i] = 3.0
        #	elif kind[i] == 'Bi': # Bi5+ 
        #		polyvalent = 1
        #		R0[i] = 2.050
        #		B[i] = 0.318
        #		V[i] = 5.0
	        elif kind[i] == 'Th': 
		        R0[i] = 2.117
		        B[i] = 0.420
		        V[i] = 4.0
        #	elif kind[i] == 'U': # U4+ 
        #		polyvalent = 1
        #		R0[i] = 2.100
        #		B[i] = 0.373
        #		V[i] = 4.0
        #	elif kind[i] == 'U': # U5+ 
        #		polyvalent = 1
        #		R0[i] = 2.009
        #		B[i] = 0.660
        #		V[i] = 5.0
	        elif kind[i] == 'U': # U6+ 
		        polyvalent = 1
		        R0[i] = 2.046
		        B[i] = 0.473
		        V[i] = 6.0
	        elif kind[i] == 'Np': # Np5+ 
		        polyvalent = 1
		        R0[i] = 2.036
		        B[i] = 0.411
		        V[i] = 5.0
        #	elif kind[i] == 'Np': # Np6+ 
        #		polyvalent = 1
        #		R0[i] = 2.022
        #		B[i] = 0.523
        #		V[i] = 6.0
        #	elif kind[i] == 'Np': # Np7+ 
        #		polyvalent = 1
        #		R0[i] = 2.076
        #		B[i] = 0.477
        #		V[i] = 7.0
	        elif kind[i] == 'Am': 
		        R0[i] = 2.068
		        B[i] = 0.392
		        V[i] = 3.0
	        elif kind[i] == 'Cm': 
		        R0[i] = 2.034
		        B[i] = 0.412
		        V[i] = 3.0
	        elif kind[i] == 'O':  # as O is anion
		        R0[i] = 0.0
		        B[i] = 0.0
		        V[i] = 2.0
        # not found in list
	        else:
                        print('Not found in build-in list!  kind= ',kind[i])
                        
                        

	        #if kind[i] != 'O':
		        #print(kind[i] + ': R0 = ' + str(R0[i]) + ' : B = ' + str(B[i]) + ' : Q = ' + str(V[i]))
        return R0 ,B 


"""
#1
def vb_analize_from_struct(lse,se):
        ### input ###
        #class?:LightStructureEnvironments.from_structure_environments   (from chemenv)
        valence_list=[]

        site_list=[]        
        
        for isite in range(len(lse.structure)):
                print(se.structure[isite].species_and_occu)
                site_list.append(se.structure[isite].species_and_occu)
        print('site_list',site_list)

        
        ### take vb paramater ###
        Ro_list,B_list=take_vbparemater()
        print(Ro_list)
        print(B_list)


                
        
        for isite in range(len(lse.structure)):

                center_symbol  =lse.structure[isite].species_and_occu
                center_position=lse.structure.sites[isite].coords
                neigbor_symbols=[]
                neighbor_positions=[]

                print("\n**isite=",isite) 
                print(se.structure[isite].species_and_occu,":",lse.structure.sites[isite].coords)

                isite_coords=lse.structure.sites[isite].coords

                print()

                if lse.coordination_environments[isite] is None:continue 
                for i in range(len(lse.coordination_environments[isite])):
                        print('ce_symbol  =',lse.coordination_environments[isite][i]['ce_symbol'])
                        print('ce_fraction= ',lse.coordination_environments[isite][i]['ce_fraction'])
                        cn=len(lse.neighbors_sets[isite][i].all_nbs_sites_indices)

                        _distance=[]
                        print("species:coords:index")
                        for inbs in range(len(lse.neighbors_sets[isite][i].neighb_sites_and_indices)):      #neib loop                   
                                coords=lse.neighbors_sets[isite][i].neighb_sites_and_indices[inbs]['site'].coords
                                index=lse.neighbors_sets[isite][i].neighb_sites_and_indices[inbs]['index']
                                coords_checkindex=lse.structure.sites[index].coords
                                _distance.append(np.linalg.norm(isite_coords - coords))
                                species=lse.structure.sites[index]._species
                                
                                print(species,":",coords,":",index)
                                
                                neigbor_symbols.append(species)
                                neighbor_positions.append(coords)

                        if not set(se.neighbors_sets[isite][cn][0].distances)==set(_distance):
                                print("!!!atention to nb coordination!!!")
                                print("        ",se.neighbors_sets[isite][cn][0].distances)
                                print("        ",_distance)

                        print()
                

                ### which factional or non fractional is this site 
                if re.findall(r'\b[A-Z][a-z]*(\d\D)*1\b',center_symbol):
                        val=vb_analize_from_site(center_symbol,center_position,neigbor_symbols,neighbor_positions)
                        valence_list.appnd(val)
                else :
                        print(re.findall(r'\b[A-Z][a-z]*(\d\D)*\d\b',center_symbol))
                        sys.exit()
                        #val=vb_analize_from_site_fractional(center_symbol,center_position,neigbor_symbols,neighbor_positions)
                        #valence_list.appnd(val)

        return valence_list

"""
