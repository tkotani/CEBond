import os ,re,shutil

def remove(pattern,out_path):
        try:os.mkdir(out_path)
        except:pass
        dlist=os.listdir()
        grep=[]
        for item in dlist:
                if re.findall(pattern,item):
                        grep.append(item)
                        shutil.move("./"+item,out_path)

        print('remove files ',grep)
#out_path='./contena/VBA_dat'
#pattern=r'\d+_out\.dat'                      
#remove(pattern,out_path)
