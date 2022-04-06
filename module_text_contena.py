import os
def text_output(text,filename,directry=None):
        cdir=os.getcwd()

        if directry==None:
                directry="contena"

        try:
                os.makedirs(cdir+"/"+directry+"/")
        except:pass

        with open(cdir+"/"+directry+"/"+filename,'w') as f:
              f.write(text)  


def text_input(filepass):
        with open(filepass) as f:
                data=f.read()
        return data
                

#ex
#text_output('ccc','ccccc.dat','test2/test2')
#print(text_input("./test2/test2/ccccc.dat"))
