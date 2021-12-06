import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
rcParams['font.size'] = 8
import matplotlib.pyplot as plt
import scipy
from scipy.stats import chisquare

def combine_sfds(df1, df2, diameter, sample1=None, sample2=None):
    if sample1 is not None:
        df1 = df1.sample(sample1, replace=True)
    if sample2 is not None:
        df2 = df2.sample(sample2, replace=True)
    df_1 = df1[df1['Effective Diameter [m]'] < diameter]
    df_2 = df2[df2['Effective Diameter [m]'] >= diameter]
    return pd.concat([df_1, df_2])

def combine_sfds_pro(dfs, diameters, sample=None):
    
    print(diameters, len(diameters))
    
    # Biggest diameter SFD sets sample size
    diam = diameters[-1]
    df = dfs[-1]
    if sample is not None:
        dfl = df[df['Effective Diameter [m]'] >= diam].sample(sample, replace=True)
    else:
        dfl = df[df['Effective Diameter [m]'] >= diam]
    
    # Now consider the other SFDs, sample
    for i in range(len(dfs)-1):
        df = dfs[len(dfs)-2-i]
        diam = diameters[len(diameters)-1-i]
        ntot = len(df)
        nlarger_in = len(df[df['Effective Diameter [m]'] >= diam])
        nlarger_out = len(dfl[dfl['Effective Diameter [m]'] >= diam])
        nsample = int(ntot*nlarger_out/nlarger_in)      # Number of samples of input required
        print('Resampling ',i,' from ',ntot,' to ',nsample)
        dfn = df.sample(nsample, replace=True)          # Re-sample the input SFD to match the output
        dfn = dfn[dfn['Effective Diameter [m]'] < diam] # Throw away the large craters
        dfl = pd.concat([dfl, dfn])                     # Combine the SFDs
        
    return dfl

def time_on_mars(minmass, nimpacts, ME_ratio=2.03):
    '''Returns the time in years that one would expect the specified number 
    of impactors striking Mars greater than the specified mass.  Based on observed
    impact rates on Earth, scaled to Mars. Uses Mars/Earth ratio of 2.03 by default.
    '''
    N_Earth = 10**(-0.926*np.log10(minmass) + 4.739)
    N_density_Earth = N_Earth/510.1E6
    N_density_Mars = ME_ratio*N_density_Earth
    N_Mars = N_density_Mars * 144.8E6
    return nimpacts/N_Mars

def n_mars_time(minmass, time, ME_ratio=2.03):
    '''Returns the number of impactors striking Mars greater than
    the specified mass in the specified time in years.  Based on observed
    impact rates on Earth, scaled to Mars. Uses Mars/Earth ratio of 2.03.
    '''
    N_3E = 10**(-0.926*np.log10(3.) + 4.739) # Number of impacts >3 kg on Earth per year
    N_3M = N_3E*(144.8/510.1)*ME_ratio
    return int(N_3M * (minmass/3.)**-0.926 * time)

def combine_sfds_mass(dfs, diameters, masses, time=1.):
    '''Function to merge separate crater SFDs based on minimum mass
        of impactors for each SFD and time in years. Merge takes
        portion of SFD divided by bounding diameters.
    '''
    assert len(dfs) == len(diameters)
    assert len(dfs) == len(masses)
    assert len(masses) == len(diameters)
    
    # Sample from input SFD based on min. mass and time
    df = dfs[-1].sample(n_mars_time(masses[-1], time), replace=True)

    # For the biggest mass, keep the large craters
    dfo = df[df['Effective Diameter [m]'] >= diameters[-1]] 

    # For the smaller masses, add the small craters
    for i in range(len(dfs)-1):
        df = dfs[len(dfs)-2-i].sample(n_mars_time(masses[len(dfs)-2-i], time), replace=True)
        dfn = df[df['Effective Diameter [m]'] < diameters[len(diameters)-1-i]] 
        dfn = dfn[dfn['Effective Diameter [m]'] >= diameters[len(diameters)-2-i]] 
        dfo = pd.concat([dfo, dfn])
          
    return dfo

def classify_clusters(df):
    '''Classify results of simulation into clusters, singulars and airbursts 
    based on the number of craters produced.'''
    clusters=df[df['No. of Craters']>1]
    singulars=df[df['No. of Craters']==1]
    bursts=df[df['No. of Craters']==0]
    return clusters, singulars, bursts

def comp_plot(model, obs, var1, var2, xlimit, ylimit=None):
    sns.color_palette("light:b", as_cmap=True)
    fig = plt.figure(figsize=(9,3))
    fig.add_subplot(132)
    plt.title('Model')
    sns.kdeplot(x=model[var1], y=model[var2], shade=True, cmap='light:b')
    plt.xlim(xlimit)
    if ylimit:
        plt.ylim(ylimit)
    fig.add_subplot(131)
    plt.title('Observations')
    sns.kdeplot(x=obs[var1], y=obs[var2], shade=True, cmap='light:b')
    plt.xlim(xlimit)
    if ylimit:
        plt.ylim(ylimit)
    fig.add_subplot(133)
    plt.title('Comparison')
    sns.kdeplot(x=obs[var1], y=obs[var2], shade=True, cmap='light:b')
    sns.kdeplot(x=model[var1], y=model[var2], shade=False)
    plt.xlim(xlimit)
    if ylimit:
        plt.ylim(ylimit)
    fig.tight_layout()
    
def comp_hist(df, var, bin_list, stat='probability', **kwargs):
    sns.histplot(df, x=var, multiple='dodge', palette='Paired',
                 element='bars', hue='Type', kde=True, stat=stat,
                 bins=bin_list, shrink=0.8, **kwargs)
    return

def comp_hist_mpl(df, var, bin_list=15, **kwargs):
    dfm = df[df['Type'] == 'Model']
    dfo = df[df['Type'] == 'Observed']
    plt.hist((dfm[var], dfo[var]), bins=np.linspace(0,1.,11), density=True,
                 color=sns.color_palette("Paired")[0:2], label=['Model','Observed'], edgecolor='black', alpha=0.6)
    sns.kdeplot(data=df, x=var, multiple='layer', hue='Type', palette='Paired', common_norm=False, clip=[0.,1.]) #, **kwargs)

def cluster_stats(oc, os):
    n_oc     = len(oc)
    n_oc_g5  = len(oc[oc['Effective Diameter [m]'] > 5.65])
    n_oc_g10 = len(oc[oc['Effective Diameter [m]'] > 10. ])
    oc = oc[oc['No. of Craters'] > 5]
    n_oc5 = len(oc)
    n_oc5_g5 = len(oc[oc['Effective Diameter [m]'] > 5.65])
    n_oc5_g10 = len(oc[oc['Effective Diameter [m]'] > 10. ])

    n_os     = len(os)
    n_os_g5  = len(os[os['Effective Diameter [m]'] > 5.65])
    n_os_g10 = len(os[os['Effective Diameter [m]'] > 10. ])

    print('                       All  : ', n_oc, n_os, round(n_oc/n_os, 2))
    print('               D_eff > 5.65 : ', n_oc_g5, n_os_g5, round(n_oc_g5/n_os_g5, 2))
    print('Nc>5 clusters, D_eff > 5.65 : ', n_oc5_g5, n_os_g5, round(n_oc5_g5/n_os_g5, 2))
    print('               D_eff > 10   : ', n_oc_g10, n_os_g10, round(n_oc_g10/n_os_g10, 2))
    print('Nc>5 clusters, D_eff > 10   : ', n_oc5_g10, n_os_g10, round(n_oc5_g10/n_os_g10, 2))
    
    return n_oc, n_oc_g5, n_oc_g10, n_oc5, n_oc5_g5, n_oc5_g10, n_os, n_os_g5, n_os_g10

def chi_square_comparison(dist_mod, dist_obs, var, nbins=10, logbins=False, **kwargs):
    ''' Compares model and observation distributions using the chi-square statistic '''

    if logbins:
        od = np.log10(dist_obs[var])
        sd = np.log10(dist_mod[var])
    else:
        od = dist_obs[var]
        sd = dist_mod[var]

    # Bin the observed and synthetic data consistently. . .
    h_obs, bins = np.histogram(od, bins=nbins, **kwargs)
    h_mod, bins = np.histogram(sd, bins=bins, **kwargs)
        
    return chisquare(f_obs=h_mod, f_exp=h_obs)

def compare_histograms_clusters(combined, ID, **kwargs):

    #combined['Large Crater Fraction'] = combined['Large Crater Fraction'] * 100
    #combined['Aspect Ratio'] = combined['Aspect Ratio'] * 100
    #combined['Aspect Ratio (alt)'] = combined['Aspect Ratio (alt)'] * 100

    fig = plt.figure(figsize=(10,8))
    fig.add_subplot(221)
    var = 'Dispersion [m]'
    comp_hist(combined, var, bin_list=15, log_scale=True, **kwargs)
    plt.xlabel(var)
    #plt.title(var+': model vs obs.')

    fig.add_subplot(222)
    var = 'Aspect Ratio (alt)'
    comp_hist_mpl(combined, var, bin_list=np.linspace(0.,1.,11), log_scale=False, **kwargs)
    plt.xlabel('Aspect Ratio [%]')
    #plt.title(var+': model vs obs.')

    fig.add_subplot(223)
    var = 'Large Crater Fraction'
    comp_hist(combined, var, bin_list=np.linspace(0.,1.,11), log_scale=False, **kwargs)
    plt.xlabel(var+' [%]')
    #plt.title(var+': model vs obs.')

    var = 'CSFD exponent'
    fig.add_subplot(224)
    comp_hist(combined, var, bin_list=10, **kwargs)
    plt.xlabel(var)
    #plt.title(var+': model vs obs.')

    plt.tight_layout()
    plt.savefig('Histograms_Clusters_'+ID+'.pdf')

def compare_histogram_SCN(combined_all, combined_lar, singles, ID, logy=False, **kwargs):

    fig = plt.figure(figsize=(15,4))
    fig.add_subplot(131)
    var = 'Effective Diameter [m]'
    comp_hist(singles, var, bin_list=np.linspace(0.,1.8062,13), log_scale=True, **kwargs)
    plt.xlabel(var)
    plt.title('Single craters')
    if logy:
        plt.yscale('log')
        plt.ylim(5.E-4,5.E-1)

    fig.add_subplot(132)
    var = 'Effective Diameter [m]'
    comp_hist(combined_all, var, bin_list=np.linspace(0.,1.806,13), log_scale=True, **kwargs)
    plt.xlabel(var)
    plt.title('Clusters')
    if logy:
        plt.yscale('log')
        plt.ylim(5.E-4,5.E-1)
        
    fig.add_subplot(133)
    var = 'No. of Craters'
    comp_hist(combined_lar, var, bin_list=15, log_scale=True, **kwargs)
    plt.xlabel(var)
    plt.title(var)
    if logy:
        plt.yscale('log')
        plt.ylim(5.E-4,5.E-1)

    plt.tight_layout()
    plt.savefig('Histograms_SCN_'+ID+'.pdf')

#### Comparison between Model and Observation
def score_model(ID=None, synthetics=pd.DataFrame(), 
                Deff_c=5.65, N_c=5, size=10, nbins=[15, 10, 15, 15, 10, 15, 15, 10], verbose=False):

    stats = {}

    # Get column headings from Obs. file
    observedclusters=pd.read_csv('./obs-data/observedclusters.csv')
    observedsingulars=pd.read_csv('./obs-data/observedsingulars.csv')
    observed = pd.concat([observedclusters, observedsingulars], sort=False)
    characteristics = observedclusters.columns.tolist()
    characteristics.remove('Minimum Diameter [m]')
    characteristics.remove('Median Diameter [m]')
    if verbose:
        print(characteristics)

    if ID is not None:
        synthetics=pd.read_csv('./model-data/output-'+ID+'.csv') #.sample(2092, replace=True)
    syntheticclusters, syntheticsingulars, syntheticbursts = classify_clusters(synthetics)
    if verbose:
        print()
        print('Observed: ')
        cluster_stats(observedclusters, observedsingulars)
        print()
        print('Model: ')
        cluster_stats(syntheticclusters, syntheticsingulars)
        print()

    # Filter the results and data to exclude small clusters and those with few craters
    syntheticclusters = syntheticclusters[syntheticclusters['No. of Craters'] > N_c]
    syntheticclusters = syntheticclusters[syntheticclusters['Effective Diameter [m]'] > Deff_c]
    observedclusters = observedclusters[observedclusters['No. of Craters'] > N_c]
    observedclusters = observedclusters[observedclusters['Effective Diameter [m]'] > Deff_c]
    
    kstests = []; estests = []; cstests = []
    for char,a in zip(characteristics[1:],range(8)):
        ocdf = observedclusters[char].value_counts().sort_index()[::-1].cumsum()/len(observedclusters[char])
        if len(syntheticclusters) > 0:
            cdf = syntheticclusters[char].value_counts().sort_index()[::-1].cumsum()/len(syntheticclusters[char])
        if (char != "CSFD exponent" and char != "Aspect Ratio" and 
            char != "Large Crater Fraction"):
            try:
                cstest, _ = chi_square_comparison(syntheticclusters, observedclusters, char,
                                                     nbins=nbins[a], logbins=True, density=True)
                kstest=scipy.stats.ks_2samp(np.log10(syntheticclusters[char]), 
                                            np.log10(observedclusters[char]),alternative='two-sided')[0]
                estest=scipy.stats.epps_singleton_2samp(np.log10(syntheticclusters[char]), 
                                            np.log10(observedclusters[char]))[0]
            except:
                kstest = 1.; estest = 100.; cstest = 100.
        else:
            try:
                cstest, _ = chi_square_comparison(syntheticclusters, observedclusters, char, 
                                                 nbins=nbins[a], logbins=False, density=True)
                kstest=scipy.stats.ks_2samp(syntheticclusters[char], 
                                        observedclusters[char],alternative='two-sided')[0]
                estest=scipy.stats.epps_singleton_2samp(syntheticclusters[char], 
                                        observedclusters[char])[0]
            except:
                kstest = 1.; estest = 100.; cstest = 100
        kstests.append(kstest); estests.append(estest); cstests.append(cstest)
        stats[char] = cstest
        if verbose:
            print('%25s %3.2f %4.1f %4.2f' % (char, kstest, estest, cstest))

    cs_ratio = (len(syntheticclusters[syntheticclusters['Effective Diameter [m]'] > size]) / 
              len(syntheticsingulars[syntheticsingulars['Effective Diameter [m]'] > size]))
    stats['C-S ratio'] = cs_ratio
    stats['Total'] = np.sum(cstests)
    if verbose:
        print('%25s %3.2f %4.3f %3.2f' % ("Total", np.sum(kstests[2:]), np.sum(estests[2:]), stats['Total']))
    
    stats['No. of Clusters'] = len(syntheticclusters)
    
    return stats
