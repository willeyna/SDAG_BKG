from imports import *

####### PARAMETERS
ntrials = 10000
ninj_t = 3
ninj_c = 5
#name of classic likelihood function to use as a baseline for comparison
LLH_name = "LLH_detector"


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


######### READING IN BACKGROUND FILES


fils=glob.glob("./outputs/*.npz")
filtered_fils = []
llh_fils = []
#filters to make sure each file to be repacked is from the right method and is finalized
for fil in fils:
    dat = np.load(fil)
    # READS IN EVERY DATA FILE IN OUTPUTS TO FIND ONES USING THE NAMED METHOD AND BKG COUNT
    if dat["method"] == Method and dat["final"] == True and (dat['Bkg'] == B).all():
        filtered_fils.append(dat["bkg_TS"])
    if dat["method"] == LLH_name and dat["final"] == True and (dat['Bkg'] == B).all():
        llh_fils.append(dat["bkg_TS"])

if len(llh_fils) == 0:
    print("ERROR: Could not find appropriate base LLH file (Try checking background values)")

TA_bkg_TS = np.sort(np.concatenate(filtered_fils).flatten(), axis = 0)
LLH_bkg_TS = np.sort(np.concatenate(llh_fils).flatten(), axis = 0)

print('Files read in')

#error if either bkg sets are empty
if len(LLH_bkg_TS) == 0 or len(TA_bkg_TS) == 0:
    print('ERROR: One of the background distributions is missing and/or not being properly read in')


######################################## SIGNAL TRIAL CREATION

###### Signal Distribution Creation #########
TA_signal_TS = np.zeros(ntrials,dtype=np.float64)
LLH_signal_TS = np.zeros(ntrials,dtype=np.float64)

for j in (np.arange(ntrials)):
    #keeps randomizing the track sky and calculating the llh at this specified point under new skies
    evs = ev_maker(ninj_t,bkg_t,ninj_c,bkg_c)

    tracks = evs[:ninj_t+bkg_t,:]
    cascades = evs[ninj_t+bkg_t:,:]


    TA_signal_TS[j] = eval(Method + '(tracks,cascades,nin_ra, nin_dec)')[0]
    LLH_signal_TS[j] = LLH_detector(tracks,cascades, nin_ra, nin_dec)[0]




######################################## TESTING/PLOTTING CODE

###### Signifigance Creation #########

#the + .001 lets me use np.max() to see overflow amount; Probably a better way to do it
LLH_sig = np.array([pd2sig(p_value(x, LLH_bkg_TS)) if p_value(x,LLH_bkg_TS) != 0 else pd2sig(1/LLH_bkg_TS.shape[0])+.001 for x in (LLH_signal_TS)])
LLH_mean, LLH_std = np.mean(LLH_sig), np.std(LLH_sig)
#gets the amount of overflow
LLH_overflow = np.sum(LLH_sig == LLH_sig.max())
#crops out the overflow events s.t. they can be plotted in a different color
plotLLH_sig = LLH_sig[:-LLH_overflow]

#creates an array of the overflow events
plotLLH_overflow = np.ones(LLH_overflow) * pd2sig(1/ LLH_bkg_TS.shape[0])

TA_sig = np.array([pd2sig(p_value(x, TA_bkg_TS)) if p_value(x, TA_bkg_TS) != 0 else pd2sig(1/ TA_bkg_TS.shape[0])+.001 for x in (TA_signal_TS)])
TA_mean, TA_std = np.mean(TA_sig), np.std(TA_sig)
TA_overflow = np.sum(TA_sig == TA_sig.max())
#crops out the overflow events s.t. they can be plotted in a different color
plotTA_sig = TA_sig[:-TA_overflow]

#creates an array of the overflow events
plotTA_overflow = np.ones(TA_overflow) * pd2sig(1/ TA_bkg_TS.shape[0])


##### Signifigance Difference Creation #####
diff = TA_sig - LLH_sig
diff_mean, diff_std = np.mean(diff), np.std(diff)



######### Plot Creation ############
'''
Adds a label in the top left corner stating mean, std, and number of overflow
'''

fig, axes = plt.subplots(1,3, figsize = (20,4))

textstr = '\n'.join((
    r'$\mu=%.2f$' % (LLH_mean, ),
    r'$\mathrm{overflow}=%.2f$' % (LLH_overflow, ),
    r'$\sigma=%.2f$' % (LLH_std, )))

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
axes[0].text(.07, .93, textstr, transform=axes[0].transAxes, fontsize=14,
        verticalalignment='top', bbox=props)



textstr = '\n'.join((
    r'$\mu=%.2f$' % (TA_mean, ),
    r'$\mathrm{overflow}=%.2f$' % (TA_overflow, ),
    r'$\sigma=%.2f$' % (TA_std, )))

axes[1].text(.07, .93, textstr, transform=axes[1].transAxes, fontsize=14,
        verticalalignment='top', bbox=props)


textstr = '\n'.join((
    r'$\mu=%.2f$' % (diff_mean, ),
    r'$\sigma=%.2f$' % (diff_std, )))

axes[2].text(1.3, .93, textstr, transform=axes[1].transAxes, fontsize=14,
        verticalalignment='top', bbox=props)

axes[0].hist(LLH_sig, alpha = 0.8, bins = np.arange(0,7,.25))
axes[0].hist(plotLLH_overflow, color = 'red', alpha = 0.8, bins = np.arange(0,7,.25))
axes[0].set_title(f'Basic LLH Signifigance ninj = [{ninj_t},{ninj_c}]; [{bkg_t}, {bkg_c}] bkg')

axes[1].hist(TA_sig, alpha = 0.8, bins = np.arange(0,7,.25))
axes[1].hist(plotTA_overflow, color = 'red', alpha = 0.5, bins = np.arange(0,7,.25))
axes[1].set_title(f'Topology-Aware Signifigance ninj = [{ninj_t},{ninj_c}]; [{bkg_t}, {bkg_c}] bkg')

axes[2].hist(diff,bins = 20)
axes[2].set_title(f'Difference (TA - LLH) ninj = [{ninj_t},{ninj_c}]; [{bkg_t}, {bkg_c}] bkg')


plt.savefig(f'./outputs/SIG_Compare_' + Method + "_" + str(datetime.datetime.now())[:10])

os.system(f"rm {Method}_VIZ.sb")
os.system(f"rm {Method}_comparison.out")
