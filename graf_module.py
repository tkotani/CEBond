
###  color_set  ###
# input single ce symbol 
# outoput single ce color 
import re ,subprocess
import numpy as np
base_colors=[(0.5,0,0)    ,(0.75,0.75,0),(1,0.5,0)    ,(0,0,0.5)       ,(0.5,0.5,0),
             (0,0.40,0)   ,(0.5,0,1)    ,(1,0.25,0.25),(0.5,0.25,0)    ,(0.5,0,0.5),
             (0,0.25,0.25),(0,0.75,0)   ,(0,0,0)      ,(0.75,0.75,0.75)             ]


## all ce to all different color
def __sym_to_color(symbol):       
        base_colors=[(0.5,0,0)    ,(0.75,0.75,0),(1,0.5,0)    ,(0,0,0.5)       ,(0.5,0.5,0),
                     (0,0.40,0)   ,(0.5,0,1)    ,(1,0.25,0.25),(0.5,0.25,0)    ,(0.5,0,0.5),
                     (0,0.25,0.25),(0,0.75,0)   ,(0,0,0)      ,(0.75,0.75,0.75)             ]

        if re.findall(r'\A\d+\Z',symbol):
                color=base_colors[int(symbol)-1]
                return color

        contents=subprocess.getoutput('grep -w {symbol} cn_symbols'.format(symbol=symbol))
        contents=contents.split()

        num1=int(contents[0])
        num2=int(contents[1])
        num3=int(contents[2])
        
        color=[0,0,0] 
        for i in range(3):
                if base_colors[num1-1][i]==1:                      
                        color[i]=base_colors[num1-1][i] -0.40  +0.40*num2/num3
                elif base_colors[num1-1][i]==0:
                        color[i]=base_colors[num1-1][i] + 0.40*num2/num3
                else:
                        color[i]=base_colors[num1-1][i] -0.20 + 0.40*num2/num3
        return color

## ver.2 adapt option more than 7 cn is base color
def _sym_to_color(symbol,cn7_base=None):       
        base_colors=[(0.5,0,0)    ,(0.75,0.75,0),(1,0.5,0)    ,(0,0,0.6)       ,(0.5,0.5,0),
                     (0,0.40,0)   ,(0.5,0,1)    ,(1,0.25,0.25),(0.5,0.25,0)    ,(0.5,0,0.5),
                     (0,0.25,0.25),(0,0.75,0)   ,(0,0,0)      ,(0.75,0.75,0.75)             ]

        if re.findall(r'\A\d+\Z',symbol):
                color=base_colors[int(symbol)-1]
                return color

        contents=subprocess.getoutput('grep -w {symbol} cn_symbols'.format(symbol=symbol))
        contents=contents.split()

        num1=int(contents[0])
        num2=int(contents[1])
        num3=int(contents[2])
        
        color=[0,0,0] 
        for i in range(3):

                ### more than 7 cn 
                if cn7_base:
                        if num1>=7:
                                color[i]=base_colors[num1-1][i]
                                continue

                if base_colors[num1-1][i]==1:                      
                        color[i]=base_colors[num1-1][i] -0.40  +0.40*num2/num3
                elif base_colors[num1-1][i]==0:
                        color[i]=base_colors[num1-1][i] + 0.40*num2/num3
                else:
                        color[i]=base_colors[num1-1][i] -0.20 + 0.40*num2/num3
        return color






# Ver.3 2020 04 13
def sym_to_color(symbol,cn7_base=None):    

        print()
        print("sym_to_color")   

        base_colors=[(0.2,0.2,0.2),(0.75,0.75,0),(0.9,0,0.2),(0,0,1),(0.6,0.4,0.1),(0,0.4,0),
                     (0.3,0.6,1),(1,0.25,0.25),(0.4,0.5,1),(0.8,0.6,1),
                     (0,0.45,0.7),(0.5,1,0.3),(0,0,0),(0.5,0.5,0.5) ]

        if re.findall(r'\A\d+\Z',symbol):  #"1", "7", "13" >> rgb
                color=base_colors[int(symbol)-1]
                return color

        contents=subprocess.getoutput('grep -w {symbol} cn_symbols'.format(symbol=symbol))
        contents=contents.split()

    
        num1=int(contents[0])# cn
        num2=int(contents[1])# ce
        num3=int(contents[2])# a number of ce
        
        base_color=list(base_colors[num1-1])
        print("basecolor=",base_color)
        print(num1,num2,num3)



        over_dic={}
        color=[0]*3
        for i in range(3):
                print("i=",i)
                if cn7_base:   ### more than 7 cn 
                        if num1>=7:
                                base_color[i]=base_colors[num1-1][i]
                                continue

                
                #if base_colors[num1-1][i]==1:                      
                #        color[i]=base_colors[num1-1][i] -0.40  +0.40*num2/num3
                #elif base_colors[num1-1][i]==0:
                #        color[i]=base_colors[num1-1][i] + 0.40*num2/num3
                #else:
                #        color[i]=base_colors[num1-1][i] -0.20 + 0.40*num2/num3
                
                 
                #color_num=base_colors[num1-1][i] -0.40+ 0.80*num2/num3



                gradation=num3*0.15

                #lower    =0.5*gradation
                #print(base_color,sum(base_color))
                #if sum(base_color)<=1:lower=0.1*gradation
                #if sum(base_color)>2:lower=0.9*gradation

                lower=sum(base_color)/3 *gradation


                print("gradation",gradation,"lower",lower)



                if num3-1==0:
                        color_num=base_colors[num1-1][i]
                else:        
                        color_num=base_colors[num1-1][i] -lower+ gradation*(num2-1)/(num3-1)
                print("cnum",color_num)



                if color_num<0:
                        over=color_num
                        print("over+",over)
                        over_dic[i]=over
                        color_num=0

                if color_num>1:
                        over=color_num-1
                        print("over-",over)
                        over_dic[i]=over
                        color_num=1

                color[i]=color_num
        print("ce_color_1=",color)
        if over_dic:
                print("over=",over_dic)
                over_rgb=over_dic.keys()

                
                safe_rgb=set([0,1,2])-set(over_rgb)
                        
                
                print("o",over_rgb)
                print("s",safe_rgb)
                
                all_over=sum([over_dic[key] for key in over_rgb])
                print(all_over)
                
                for safe in safe_rgb:
                        color[safe]+=all_over/len(safe_rgb)
                
        for i in range(3):
             if color[i]<0:color[i]=0   
             if color[i]>1:color[i]=1

        """
        i=num2%3
        color_num=base_colors[num1-1][i] -0.40+ 0.80*(num2-1)/(num3-1)

        if color_num<=0:
                color_num=0
        if color_num>=1:
                color_num=1
        color[i]=color_num
        """
        print(symbol,color)
      
        return color


def inversion(color):
        new_color=[0,0,0]
        #for i in range(3):
        #        new_color[i]=1-color[i]
        
        if 1.>sum(color):new_color=[0.5,0.5,0.5]
        else:   new_color=[0,0,0]  
   
        print(color,">>>>",new_color)
        return new_color



### distance Discretization
def distanse_divider(distance,class_width):
        n=class_width
        i=0
        while n<distance:
                n=n+class_width
                i=i+1
        return i


def distanse_divider_log(distance,class_width):

        log_dis=np.log(distance)

        n=class_width
        i=0
        while n<log_dis:
                n=n+class_width
                i=i+1


        return i


## sort ce label
##2019 12 24
import sys
def _sort_ce(labels):
        # sort coordinaition numbar
        sorted_labels=[]
        for cn in range(1,14):
                for ce in labels: 
                        if re.findall('.*:{d}\Z'.format(d=cn),ce):
                                sorted_labels.append(ce)
                                
                        elif ce==str(cn) :
                                sorted_labels.append(ce)

        if not len(labels)==len(sorted_labels):
                print("!!! failed sort label")
                sys.exit()
        return sorted_labels


# Ver.2 import cn_symbols
def sort_ce(labels):
        with open('cn_symbols') as l:
                cn_symbols=l.read().split('\n')

        _model_ce_list=[]
        for line in cn_symbols:
                if line=="":continue
                _model_ce_list.append(line.split(" ")[3])

        model_ce_list=[]
        cn=1
        for ce in _model_ce_list:
                if ":" in ce:
                        current_cn=ce.split(":")[1]
                else:current_cn=ce

                if cn==current_cn:
                       cn=current_cn
                       model_ce_list.append(ce)
                else:
                       model_ce_list.append(cn)
                       cn=current_cn
                       model_ce_list.append(ce)        

        # sort coordinaition numbar
        sorted_labels=[]
        for symbol in model_ce_list:
                for ce in labels: 
                        if symbol==str(ce):
                                sorted_labels.append(ce)
                        elif ce==str(cn) :
                                sorted_labels.append(ce)

        if not len(labels)==len(sorted_labels):
                print("!!! failed sort label")
                print(len(labels),len(sorted_labels),"\n",labels,"\n",sorted_labels)
                sys.exit()
        return sorted_labels


"""
def 
        for cn in range(1,14):
                for key in contena[cation].keys(): 
                        re_fd='.*:{d}\Z'.format(d=cn)
                        if re.findall(re_fd,key):
                                label.append(key)
                        if key==str(cn) :
                                label.append(key)
        print("label= ",label)
"""



import numpy as np
def light_grad_to_rgb(hue,lightness,intensity):

        #hue purple=0 [1,0,1] to red purple=1
        #lightness=0~1
        #intensity=0~1


        ### hue
        #sum(rgb_hue)=1
        purple_hue=np.array([0.5,0,0.5]) #not rgb
        rgb_hue=purple_hue

        phase=0.01/3
        
        """
        rb_role=np.array([[np.cos(2*np.pi*hue),0,-np.sin(2*np.pi*hue)],
                          [0,0,0],
                          [np.sin(2*np.pi*hue),0,np.cos(2*np.pi*hue)]])

        rgb_hue=np.dot(rb_role,purple_hue)
        print("hue=",rgb_hue)
        """
        i=0+phase
        print("#1")
        while i<1/6 and i<hue:
                rgb_hue[0]+=-phase*6*0.5
                rgb_hue[2]+=phase*6*0.5
                i+=phase
                #print(rgb_hue)

        rgb_hue=[0,0,1]
        print("#2")
        while i<3/6 and i<hue:
                rgb_hue[2]+=-phase*6*0.5
                rgb_hue[1]+=phase*6*0.5
                i+=phase
                #print(rgb_hue)

        rgb_hue=[0,1,0]
        print("#3")
        while i<5/6 and i<hue:
                rgb_hue[1]+=-phase*6*0.5
                rgb_hue[0]+=phase*6*0.5
                i+=phase
                #print(rgb_hue)
        
        rgb_hue=[1,0,0]
        print("#4")
        while i<6/6 and i<hue:
                rgb_hue[0]+=-phase*6*0.5
                rgb_hue[2]+=phase*6*0.5
                i+=phase
                #print(rgb_hue)
        print("hue=",rgb_hue)
        #rgb=round(rgb,3)

        ### lightness
        rgb_lightness=[0]*3

        for i in range(3):
                rgb_lightness[i]=rgb_hue[i]*lightness/0.333
        print("light1=",rgb_lightness)
        for i in range(3):
                if rgb_lightness[i]>1:
                        rgb_lightness[(i-1)%3]+=(rgb_lightness[i]-1)/2
                        rgb_lightness[(i+1)%3]+=(rgb_lightness[i]-1)/2
                        rgb_lightness[i]=1
                
        print("lightness2=",rgb_lightness)
        #srgb_hue

        

        #return rgb
print(light_grad_to_rgb(0.4,1,1))  



#def shannon_arrow


