imports.py 
Imported by bkg_maker.py; includes any created functions to be used 
Event sampling and LLH functions to test
Do not put anything slow in here 
Any LLH function must return TS as return[0] (and obviously as float)
LLH functions must take LLH(tracks, cascades, in_ra = default, in_dec = default)
Events of shape [nev, 5] [ra, dec, ang_error, Energy, Topology (0=track, 1 = cascade)]

create.py
Creates and runs a slurm DAG (directed acyclic graph) which manages the creation of your bkg data
Creates repack.sb with proper information 
Feed in (“LLH_name”, N_trials, bkg_num [ntrack_bkg, ncasc_bkg], Njob)
Need sdag.py in directory and sdag in PATH

bkg_maker.py
Gets submitted by jobs to create mini background npz
Creates many small npz files with N/Njob trials in each

repackage.py
Reads all npz files in cd and checks for method and final == False
Stacks these npz into one 
Removes old npz files, .out() files, submission files, and manager file 
