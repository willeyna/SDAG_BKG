from imports import *

#method function name 
Method = str(sys.argv[1])
#ntrials per job
ntrials = int(sys.argv[2])
#background number
B = np.array(sys.argv[3].split(','), dtype = 'int')
#job number
tag = int(sys.argv[4])
#total bkg trials
N = int(sys.argv[5])

ntrack, ncasc = B[0], B[1]

bkg_TS = np.zeros(ntrials,dtype=np.float64)

t1= datetime.datetime.now()
for j in (np.arange(ntrials)):
    #creates track and cascade background events
    ev = ev_maker(0,ntrack,0,ncasc)
    tracks = ev[:ntrack]
    cascades = ev[ntrack:]

    bkg_TS[j] = eval(Method + '(tracks,cascades)')[0]

t2= datetime.datetime.now()
dt=(t2-t1).total_seconds()

#saves the bkg sample
#method allows the method to be checked in repacking and final tells repackager whether to repack it (if final is False) or not (True)
np.savez(f'{Method}_bkg_{tag}', bkg_TS = bkg_TS, method = Method, date= str(datetime.datetime.today())[:10], dt=dt, final = False)
