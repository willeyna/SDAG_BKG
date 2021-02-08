import glob
import numpy as np
import socket
import datetime
import os
import sys

Method = str(sys.argv[1])


fils=glob.glob("./*.npz")
filtered_fils = []
#filters to make sure each file to be repacked is from the right method and is not finalized
for fil in fils:
    dat = np.load(fil)
    if dat["method"] == Method and dat["final"] == False:
        filtered_fils.append(fil)
        
fils = filtered_fils
#filters 
todo=len(fils)
print(f"Files found: {todo}")
outdat=np.load(fils[0])


#baby date string
datestr = datetime.datetime.now()
bb_datestr = f"{datestr.month}-{datestr.day}-{datestr.year}"

outdat=outdat["bkg_TS"]
os.system('rm ' + fils[0])
os.system('rm '+ Method + 'bkg0' + '.out')
os.system('rm run_bkg0' + '.sb')



statstr=" "
if(todo):
    for j,fil in enumerate(fils[1:]):
        print(len(statstr)*" ",end="\r")
        statstr=f"merging file {j+1}, {fil}"
        print(statstr,end="\r")
        file=np.load(fil)
        outdat=np.hstack((outdat,file["bkg_TS"]))
        os.system('rm ' + fil)
        os.system('rm '+ Method + 'bkg' + str(j+1) + '.out')
        os.system('rm run_bkg' + str(j+1) + '.sb')

os.system('rm ' + Method + '_manager.sdag')
os.system('rm ' + Method + '_repack.sb')

print("sorting data")

outdat.sort(axis=0)

print("packing data")

#rename
namestr="LLH_STACK_bkg_TS"


filnam='./outputs/' + namestr+"_"+bb_datestr
np.savez(filnam,bkg_TS = outdat, method = Method, final = True)
print(f"file saved to {filnam}.npz")
