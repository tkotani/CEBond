#!/usr/bin/env python3
import sys,os
f= sys.argv[1]
import urllib.request

ciflist=open(f)
for i,icif in enumerate(ciflist):
    if i==0: continue
    aaa=icif.split('\t')[0]+'.cif'
    url = 'http://www.crystallography.net/cod/'+aaa
    print(url)
    if os.path.isfile(aaa) :
        continue
    else:
        with urllib.request.urlopen(url) as u:
            with open(aaa, 'bw') as ocif:
                ocif.write(u.read())
