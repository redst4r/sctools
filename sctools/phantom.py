import pandas as pd
import numpy as np
from scipy.stats import binom
import matplotlib.pyplot as plt
import tqdm
import toolz
import pandas as pd
from scipy.special import logsumexp
import scipy.special as sps
import itertools
RUSTPHANTOM_FOLDER = '/home/mstrasse/TB4/rust_code/rustphantompurger'
RUST_ENV = 'RUSTUP_HOME=/home/mstrasse/TB4/rust_installation/.rustup CARGO_HOME=/home/mstrasse/TB4/rust_installation/.cargo '

def named_files_to_argstring(filedict):
    return " ".join([f"{k}:{v}" for k,v in filedict.items()])

class PhantomRustWrapper():
    def __init__(self, rust_phantom_folder, rust_env):
        self.rustfolder=rust_phantom_folder
        self.ENV=rust_env

    def rust_phantom_wrapper(self, outfile, folders:dict, t2g):
        folderstring = named_files_to_argstring(folders)
        cmd = f'''cd {self.rustfolder}; {self.ENV} cargo run --release -- \\
-o {outfile} phantom-fingerprint --t2g {t2g} \\
--infolders "{folderstring}"'''
        return cmd

    def rust_phantom_detect_cell_overlap_wrapper(self, outfile, folders):
        folderstring = named_files_to_argstring(folders)
        cmd = f'''cd {self.rustfolder}; {self.ENV} cargo run --release -- \\
-o {outfile} phantom-cb \\
--infolders "{folderstring}"'''
        return cmd

    def rust_phantom_filter_wrapper(self,infolders:dict, outfiles:dict, outfiles_removed:dict, outfiles_ambiguous:dict, csv:str, t2g:str):

        cmd =  f"""cd {self.rustfolder}; {self.ENV} cargo run --release -- \\
--output /dev/null phantom-filter --t2g {t2g} \\
--infolders "{named_files_to_argstring(infolders)}" \\
--outfiles "{named_files_to_argstring(outfiles)}" \\
--removed "{named_files_to_argstring(outfiles_removed)}" \\
--ambiguous "{named_files_to_argstring(outfiles_ambiguous)}" \\
--csv {csv} --threshold 0.99
        """
        return cmd

def parse_cb_overlap(df_cb):
    """
    turns the output of rust phantom-cb into a dataframe nice for pairwise plotting
    i.e. it contains all pairs of samples, and the respective amount of umis of the same cell in both samples
    """
    samplenames = df_cb.drop('CB', axis=1).columns.values

    df_paired = []

    for s1, s2 in itertools.combinations(samplenames, 2):
        df_tmp = pd.DataFrame({'numi1': df_cb[s1].values, 'numi2':df_cb[s2].values}).query('numi1>0 and numi2>0')
    #     df_tmp['CB'] = cb_overlap['CB']
        df_tmp['sample1'] = s1
        df_tmp['sample2'] = s2
        df_paired.append(df_tmp)
    #     break
    df_paired = pd.concat(df_paired)
    df_paired = df_paired.query('numi1>0 and numi2>0')

    df_paired['sample1'] = pd.Categorical(df_paired['sample1'], samplenames)
    df_paired['sample2'] = pd.Categorical(df_paired['sample2'], samplenames)

    return df_paired

class Phantom():
    def __init__(self, df_phantom):
        self.df = df_phantom.copy()
        self.samplenames = df_phantom.drop('frequency', axis=1).columns.values
        self.S = len(self.samplenames)
        R=self.df[self.samplenames].sum(1) # amplifications across experiments
        K=(self.df[self.samplenames]>0).sum(1) # n-samples, if==1 its a nonchimera
        self.df['r'] = R
        self.df['k'] = K

        # the molecular complexity profile
        self.vrs = None

    def get_df_chimer(self):
         # a table of r, total molecules and non-chimer molecules
        # needed for paramtere est
        df_chimer = pd.crosstab(self.df.r, self.df.k==1, values=self.df.frequency, aggfunc=np.sum).replace(np.nan, 0)
        df_chimer['total'] = df_chimer.sum(1)
        if True not in df_chimer.columns.values:
            df_chimer[True] = 0
        if False not in df_chimer.columns.values:
            df_chimer[False] = 0

        df_chimer['n_nonchimeric']= df_chimer[True]
        df_chimer.drop([False, True], axis=1, inplace=True)
        df_chimer.reset_index(inplace=True)

        return df_chimer

    def estimate_sihr(self, rmin=1, rmax=30, plotting=False):
        df_chimer = self.get_df_chimer()
        r = df_chimer.query('r>=@rmin and r<=@rmax').r
        mr = df_chimer.query('r>=@rmin and r<=@rmax').total.values
        zr = df_chimer.query('r>=@rmin and r<=@rmax').n_nonchimeric.values
        df_binom_fit =binomial_model_fit(r, mr, zr)
        df_binom_fit_exact =binomial_model_fit_exact(r, mr, zr, self.S)
        df_binom_fit_exact = df_binom_fit_exact.sort_values('logp', ascending=False)
        phat = df_binom_fit_exact.p.values[0]

        if plotting:
            plt.figure()
            plt.scatter(np.log10(df_binom_fit.p), df_binom_fit.logp, label='approx')
            plt.scatter(np.log10(df_binom_fit_exact.p), df_binom_fit_exact.logp, s=1, label='exact')
            plt.ylim([-1e8, 12])
            plt.legend()

            plt.figure()
            plt.scatter(df_chimer.r, df_chimer['n_nonchimeric']/df_chimer['total'])
            plt.xlim([0, rmax])
            plt.ylim([0.5, 1])
            plt.plot(df_chimer.r, (1-phat)**df_chimer.r)

        return phat, df_binom_fit_exact

    def get_vr_hat(self, samplename):
        # the fraction of READ count in samplename as a function of r
        # returns a series indexed by r
        #
        # in the paper they call it molecular complexity profile
        # for fixed r, summed over samples == 1

        vr = []
        r_range = np.sort(self.df.r.unique())
#         r_range = range(1,100)

        for r in r_range:
            df_tmp = self.df.query('r==@r')
            # taking into account the frqe of observing this fingerprint
            df_tmp =df_tmp[self.samplenames].values * df_tmp['frequency'].values.reshape(-1,1)
            df_marginal = df_tmp.sum(0)
            df_marginal = df_marginal / df_marginal.sum()

            vr_value =  df_marginal[self.samplenames==samplename]
            assert len(vr_value) == 1
            vr_value =vr_value[0]
            vr.append({'r': r, 'vr': vr_value})
        vr = pd.DataFrame(vr).set_index('r')
        vr.rename({'vr': samplename}, axis=1, inplace=True)
        return vr[samplename]

    def get_all_vr(self):
        v = []
        for s in tqdm.tqdm(self.samplenames):
            vrs = self.get_vr_hat(s)
            v.append(vrs)
        self.vrs = pd.DataFrame(v)

#     def get_pi_r_hat(self, samplename, p):
# #         v = self.get_vr_hat(samplename)
#         assert self.vrs is not None
#         v = self.vrs.loc[samplename]
#         pi= (v * (self.S-1) + (p-1)) / (self.S *p -1)

#         # substitution negtive values
#         pi = np.clip(pi, 1e-6, None)
#         pi = pi/pi.sum()
#         return pi
    def get_pi_r_hat(self, p_nohop):
#         v = self.get_vr_hat(samplename)
        assert self.vrs is not None
        assert p_nohop> 0.5, "p_nohop < 0.5, make sure you didnt input the SIHR instead"
        v = self.vrs
        pi= (v * (self.S-1) + (p_nohop-1)) / (self.S *p_nohop -1)

        # substitution negtive values
        pi = np.clip(pi, 1e-6, None)
        pi = pi/pi.sum(0)
        return pi


    def posterior_via_multinomial(self, fingerprint, p_nohop ):
        """
        calcualtes the posterior explicitly using the multinomial (see page 22)

        it does yield exactly the same as posterior_sample_of_origin()
        it's mostly educational
        """

        assert p_nohop> 0.5, "p_nohop < 0.5, make sure you didnt input the SIHR instead"
        r = np.sum(list(fingerprint.values()))
        pir = self.get_pi_r_hat(p_nohop)[r].to_dict() # dict with samplenames

        post_unnorm = {s: 
            multinomial_loglikelihood(
                fp=fingerprint,
                pvec= self.get_multinomial_vector(p_nohop, s)
            ) + np.log(pir[s])
            for s in self.samplenames
        }
        # norm constant
        C = logsumexp(np.array(list(post_unnorm.values())))
        logpost = toolz.valmap(lambda x: x-C, post_unnorm)

        posterior = toolz.valmap(lambda x: np.exp(x), logpost)
        posterior = pd.Series(posterior, index=self.samplenames)

        return posterior

    def get_multinomial_vector(self, p_nohop, sample_of_origin):
        """
        for the mutlinomial vector of the index hopping matrix
        one element is p_nohop, all others are (1-p)/(S-1)
        """
        S = len(self.samplenames)
        vec = {}
        for sam in self.samplenames:
            if sam == sample_of_origin:
                vec[sam] = p_nohop
            else:
                vec[sam] = (1-p_nohop)/(S-1)
        return vec

    def posterior_sample_of_origin(self, fingerprint, p_nohop):
        assert p_nohop> 0.5, "p_nohop < 0.5, make sure you didnt input the SIHR instead"

        r = np.sum(list(fingerprint.values()))
        # nom = []
        lognom = []

        pi = self.get_pi_r_hat(p_nohop)

        for s in self.samplenames:
            pi_rs = pi.loc[s, r]
            # x = pi_rs * (self.S-1)/(1/p -1)**fingerprint[s]
            # nom.append(x)
            logx = np.log(pi_rs) + fingerprint[s] * (np.log(self.S-1) - np.log(1/p_nohop-1))
            lognom.append(logx)

        lognom = np.array(lognom)
        C = logsumexp(lognom)
        lognom -= C
        lognom = pd.Series(lognom, index=self.samplenames)

        return np.exp(lognom)

    def plot_amplification_profile(self):
        if self.vrs is None:
            self.get_all_vr()
        for s in self.samplenames:
            vr = self.vrs.loc[s]

            plt.plot(vr, label=s)
        plt.legend()
        plt.xlabel('r')
        plt.xlim([0, 100])

def binomial_model_fit(r, mr, zr):
    prange = np.logspace(-6,0,1000)
    res = []
    for p in prange:
        logp = binom(n=mr, p=(1-p)**r).logpmf(zr).sum()
        res.append({'p':p, 'logp':logp})
    return pd.DataFrame(res)

def binomial_model_fit_exact(r, mr, zr, S):

    prange = np.logspace(-6,0,1000)
    res = []
    for p in prange:
        logp = binom(n=mr, p=(1-p)**r + (S-1)*((p)/(S-1))**r).logpmf(zr).sum()
        res.append({'p':p, 'logp':logp})
    return pd.DataFrame(res)

def logfactorial(n):
    #factorial and gamma are shifted by +1 and so is gammaln and logfactorial
    return sps.gammaln(n+1)

def logbinomial(N, k):
    assert N>=k
    return logfactorial(N) - logfactorial(k) - logfactorial(N-k)

def log_multinomial(k_vector):
    N = np.sum(k_vector)
    nom = logfactorial(N)
    denom = logfactorial(k_vector).sum()
    return nom-denom

def multinomial_loglikelihood(fp, pvec):
    loglike = 0
    for s in fp.keys():
        loglike+= fp[s] * np.log(pvec[s])
    loglike-= log_multinomial(np.array(list(fp.values())))
    return loglike

# testing
#df_test = pd.DataFrame({'s1': [1, 1], 's2': [0,1], 'frequency': [1,1]})

#P = Phantom(df_test)
#P.estimate_p(plotting=True)
