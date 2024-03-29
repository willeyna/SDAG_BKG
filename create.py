import os
import numpy as np
import sys
import datetime

LLH = input("Please enter the name of the LLH function you would like to use:\n")
N = int(input("Number of background trials to run:\n"))
B = input("Enter your background value ({n_track, n_casc]?):\n")
nrun = int(input("How many jobs would you like to split this into?:\n"))
islong = None; 
while True:
    long_input = input("Enter 1 for fast-track jobs and 0 for uncapped jobs:\n")
    if long_input == "1":
        islong = True
        break;
    if long_input == "0":
        islong = False
        break;

mode = input("If you would like to create a TC space for the method, enter 'TC':\nIf you would like to create a significance comparison to LLH, enter 'COMP':\n If you would like to do neither, press Enter.\n")
if mode == "":
    mode = None


nper= N//nrun
outpath="./"
taskname="bkg"
date=datetime.datetime.now()
def mkstr2(task,i,nper,fast_track):
    tasknm=LLH+task+str(i)
    if fast_track:
        bigstr=f"""#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
    
#SBATCH --time=4:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --mem=100G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name {tasknm}      # you can give your job a name for easier identification (same as -J)
#SBATCH --output /mnt/research/IceCube/willey/SW_BKG/{tasknm}.out
    
########## Command Lines to Run ##########

cd /mnt/research/IceCube/willey/SW_BKG
export PATH=$PATH:/mnt/research/IceCube/willey/conda3/bin/
conda activate base

/mnt/research/IceCube/willey/conda3/bin/python3 /mnt/research/IceCube/willey/SW_BKG/bkg_maker.py {LLH} {nper} {B} {i} {N}"""
    else: 
        bigstr=f"""#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
    
#SBATCH --time=80:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --mem=100G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name {tasknm}      # you can give your job a name for easier identification (same as -J)
#SBATCH --output /mnt/research/IceCube/willey/SW_BKG/{tasknm}.out
    
########## Command Lines to Run ##########

cd /mnt/research/IceCube/willey/SW_BKG
export PATH=$PATH:/mnt/research/IceCube/willey/conda3/bin/
conda activate base

/mnt/research/IceCube/willey/conda3/bin/python3 /mnt/research/IceCube/willey/SW_BKG/bkg_maker.py {LLH} {nper} {B} {i} {N}"""
    return bigstr
filnam="run_"+taskname
names=[]
for i in range(0, nrun):
    filstr=filnam+str(i)+".sb"
    fout=open(outpath+filstr,"w")
    fout.write(mkstr2(taskname,i,nper, islong))
    fout.close()
    names.append(filstr)

################################################################################
# REPACK.SB MAKING

fout = open(outpath + f'{LLH}_repack.sb', 'w')
sbcontent = f"""#!/bin/bash --login

########## SBATCH Lines for Resource Request ##########

#SBATCH --time=4:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --mem=100G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name {LLH}_repacking      # you can give your job a name for easier identification (same as -J)
#SBATCH --output /mnt/research/IceCube/willey/SW_BKG/{LLH}_repacking.out

########## Command Lines to Run ##########

cd /mnt/research/IceCube/willey/SW_BKG
export PATH=$PATH:/mnt/research/IceCube/willey/conda3/bin/
conda activate base

/mnt/research/IceCube/willey/conda3/bin/python3 /mnt/research/IceCube/willey/SW_BKG/repackage.py {LLH}

"""
fout.write(sbcontent)
fout.close()

################################################################################

################################################################################
# GRID_COMP.SB MAKING
if mode == "COMP":
    fout = open(outpath + f'{LLH}_VIZ.sb', 'w')
    sbcontent = f"""#!/bin/bash --login

    ########## SBATCH Lines for Resource Request ##########

    #SBATCH --time=4:00:00             # limit of wall clock time - how long the job will run (same as -t)
    #SBATCH --mem=500G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
    #SBATCH --job-name {LLH}_comparison      # you can give your job a name for easier identification (same as -J)
    #SBATCH --output /mnt/research/IceCube/willey/SW_BKG/{LLH}_comparison.out

    ########## Command Lines to Run ##########

    cd /mnt/research/IceCube/willey/SW_BKG
    export PATH=$PATH:/mnt/research/IceCube/willey/conda3/bin/
    conda activate base

    /mnt/research/IceCube/willey/conda3/bin/python3 /mnt/research/IceCube/willey/SW_BKG/grid_comp.py {LLH} {B}

    """
    fout.write(sbcontent)
    fout.close()

################################################################################

################################################################################
# GRID_COMP.SB MAKING
if mode == "TC":
    fout = open(outpath + f'{LLH}_VIZ.sb', 'w')
    sbcontent = f"""#!/bin/bash --login

    ########## SBATCH Lines for Resource Request ##########

    #SBATCH --time=4:00:00             # limit of wall clock time - how long the job will run (same as -t)
    #SBATCH --mem=500G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
    #SBATCH --job-name {LLH}_comparison      # you can give your job a name for easier identification (same as -J)
    #SBATCH --output /mnt/research/IceCube/willey/SW_BKG/{LLH}_comparison.out

    ########## Command Lines to Run ##########

    cd /mnt/research/IceCube/willey/SW_BKG
    export PATH=$PATH:/mnt/research/IceCube/willey/conda3/bin/
    conda activate base

    /mnt/research/IceCube/willey/conda3/bin/python3 /mnt/research/IceCube/willey/SW_BKG/TC_space.py {LLH} {B}

    """
    fout.write(sbcontent)
    fout.close()

################################################################################





# SDAG MAKING
bkgjobs = ''
jstr = ''

for i in range(nrun):
    bkgjobs += f'JOB J{i} run_bkg{i}.sb\n'
    jstr += f'J{i} '

#creates sdag file contents
if mode != None:
    contents = f'JOB R {LLH}_repack.sb\n' + f'JOB C {LLH}_VIZ.sb\n' + bkgjobs + '\nPARENT ' + jstr + 'CHILD R' + '\nPARENT R CHILD C'
else:
    contents = f'JOB R {LLH}_repack.sb\n' + bkgjobs + '\nPARENT ' + jstr + 'CHILD R'

fout=open(outpath+LLH+'_manager.sdag',"w")
fout.write(contents)
fout.close()


os.system(f'sdag {LLH}_manager.sdag')
