from imports import *

####### PARAMETERS

# number of trials per combination 
# increases run time on the order of size^2 
ntrials = 300

############# SETUP

#"LLH" function name
Method = str(sys.argv[1])
#BKG number
B = np.array(sys.argv[2].split(','), dtype = 'int')
#bkg track and cascade count
bkg_t, bkg_c = B[0], B[1]

in_ra,in_dec =45, 60
#gets the index for the specified point in the sky
targ_ind=hp.ang2pix(NSIDE,in_ra,in_dec,lonlat=True)
#gets the angle back in case the angle is slightly different
nin_ra,nin_dec=hp.pix2ang(NSIDE,targ_ind,lonlat=True)


#### FILE READ IN

fils=glob.glob("./outputs/*.npz")
filtered_fils = []
#filters to make sure each file to be repacked is from the right method and is finalized
for fil in fils:
    dat = np.load(fil)
    try:
        if dat["method"] == Method and dat["final"] == True:
            filtered_fils.append(dat["bkg_TS"])
    except:
        pass;

bkg = np.sort(np.concatenate(filtered_fils).flatten(), axis = 0)

#error if bkg is empty
if len(bkg) == 0:
    print('ERROR: The background distribution is missing and/or not being properly read in')
    
    
    

def TC_space(llh_func, bkg, bkg_t, bkg_c, size = 10, ntrials = 100):
    
    space = np.zeros([size,size])
    
    #TC space signifigance comparison
    for ninj_t in range(size):
        for ninj_c in range(size):
            LLH_signal_TS = np.zeros(ntrials)
            for n in range(ntrials):
                #event creator chunk 
                evs = ev_maker(*[ninj_t, bkg_t ,ninj_c, bkg_c])
                tracks = evs[:bkg_t+ninj_t, :]
                cascades = evs[bkg_t+ninj_t:, :]

                LLH_signal_TS[n] = eval(Method + '(tracks,cascades,nin_ra, nin_dec)')[0]


            sig = np.array([pd2sig(p_value(x,bkg)) if p_value(x,bkg) != 0 else pd2sig(1/bkg.shape[0])+.001 for x in (LLH_signal_TS)])

            space[ninj_t, ninj_c] = np.sum([sig >= 5])/sig.shape[0]

    return space

TC = TC_space(Method, bkg, bkg_t, bkg_c, ntrials = ntrials)

plt.contourf(TC, levels = np.linspace(0,1,11))
plt.colorbar()
plt.plot(np.linspace(0,9,10000), (1/2.4)*np.linspace(0,9,10000))
plt.title(Method + ' TC Space')
plt.savefig("./outputs/" + Method + "_TC" + str(datetime.datetime.today())[6:10] + ".png")

os.system(f"rm {Method}_VIZ.sb")
os.system(f"rm {Method}_comparison.out")
