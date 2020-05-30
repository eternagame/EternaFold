import numpy as np
import sys

dat = np.loadtxt(sys.argv[1],skiprows=1,usecols=(5,6,7,8))

se_no_lig = (dat[:,0] - dat[:,1])**2
se_w_lig = (dat[:,2] - dat[:,3])**2

def bootstrap_inds(x):
    return np.random.choice(range(len(x)),size=len(x))

bs_vals_lig=[]
bs_vals_nolig=[]
for i in range(100):
	bs_inds = bootstrap_inds(dat)
	bs_dat = dat[bs_inds]
	no_lig_corr = np.corrcoef(bs_dat[:,0],bs_dat[:,1])[0][1]
	lig_corr =  np.corrcoef(bs_dat[:,2],bs_dat[:,3])[0][1]
	bs_vals_lig.append(lig_corr)
	bs_vals_nolig.append(no_lig_corr)

#print('RMSE no lig:\t%.2f\twith lig:\t%.2f' % (np.sqrt(np.mean(se_no_lig)), np.sqrt(np.mean(se_w_lig))))
print('Corr with lig:\t%.2f\t%.2f\t no lig:\t%.2f\t%.2f' % (np.mean(bs_vals_lig), np.std(bs_vals_lig), np.mean(bs_vals_nolig), np.std(bs_vals_nolig)))
