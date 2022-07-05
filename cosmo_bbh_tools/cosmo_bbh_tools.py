import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import h5py as h5 
from scipy import interpolate
from scipy import stats
import json


# this is a little function that we will use to make the plots more beautiful (bigger ticks, labels)
# However, you do not have to use this (just uncommoment "layoutAxes" everywhere)
from matplotlib.ticker import (FormatStrFormatter,
                               AutoMinorLocator)

def layoutAxes(ax, nameX='', nameY='', \
               labelSizeMajor = 10, fontsize = 18, second=False, labelpad=None, setMinor=True):
    """
    Tiny code to do the layout for axes in matplotlib
    """
    tickLengthMajor = 10
    tickLengthMinor = 5
    tickWidthMajor  = 1.5
    tickWidthMinor  = 1.5
    
    #rc('axes', linewidth=2)
    #label1 always refers to first axis not the twin 
    if not second:
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
            #tick.label1.set_fontweight('bold')
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
            #tick.label1.set_fontweight('bold')
    if second:
        for tick in ax.xaxis.get_major_ticks():
            tick.label2.set_fontsize(fontsize)
            #tick.label1.set_fontweight('bold')
        for tick in ax.yaxis.get_major_ticks():
            tick.label2.set_fontsize(fontsize)
            #tick.label1.set_fontweight('bold')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.2)
    ax.tick_params(length=tickLengthMajor, width=tickWidthMajor, which='major')
    ax.tick_params(length=tickLengthMinor, width=tickWidthMinor, which='minor')
    ax.set_xlabel(nameX, fontsize=fontsize,labelpad=labelpad)#,fontweight='bold')
    ax.set_ylabel(nameY, fontsize=fontsize,labelpad=labelpad)#, fontweight='bold')    
    
    if setMinor==True:
        # add minor ticks:
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())

    return ax

def get_index(n, met, metallicities):
    assert met < len(metallicities), f"metallicity index cannot exceed {len(metallicities)-1}"
    return n*len(metallicities) + met

def mask_type(data, st_type1, st_type2, dco_type='bbh'):
    if dco_type=='bbh':
        st_type_sum = 28
    if dco_type=='bns':
        st_type_sum = 26
    if dco_type=='bhns':
        st_type_sum = 27
    
    masked_data = []

    for i in range(len(data)):
        if len(data[i])==0:
            masked_data.append(data[i])
        else:
            st1 = st_type1[i]
            st2 = st_type2[i]

            mask = np.array(st1 + st2==st_type_sum)
            masked_data.append(data[i][mask])
    return masked_data




def load_data(n_systems, metallicities, work_dir):    
    Z_ALL = []
    
    M1_ZAMS_ALL = []
    M2_ZAMS_ALL = []
    
    M1_CO_ALL = []
    M2_CO_ALL = []
    
    SMA_CO_ALL = []
    ECC_CO_ALL = []
    
    TYPE1_CO_ALL = []
    TYPE2_CO_ALL = []
    
    DELAY_TIMES_ALL = []
    
       
    
    for n_system in n_systems:
        for met in metallicities:
            path = work_dir + f'n_{n_system:.2e}/met_{met:.2e}_combined.h5'
            print("loading data from", path)
            
            fdata = h5.File(path, 'r')
            
            M1_ZAMS = fdata['BSE_System_Parameters']["Mass@ZAMS(1)"][()]
            M2_ZAMS = fdata['BSE_System_Parameters']["Mass@ZAMS(2)"][()]
            
            M1_CO = []
            M2_CO = []
            
            SMA_CO = []
            ECC_CO = []

            TYPE1_CO = []
            TYPE2_CO = []

            DELAY_TIMES = []
            Z = []
                
            try:
                M1_CO = fdata['BSE_Double_Compact_Objects']["Mass(1)"][()]
                M2_CO = fdata['BSE_Double_Compact_Objects']["Mass(2)"][()]              
                
                SMA_CO = fdata['BSE_Double_Compact_Objects']["SemiMajorAxis@DCO"][()]
                ECC_CO = fdata['BSE_Double_Compact_Objects']["Eccentricity@DCO"][()]
                TYPE1_CO = fdata['BSE_Double_Compact_Objects']["Stellar_Type(1)"][()]
                TYPE2_CO = fdata['BSE_Double_Compact_Objects']["Stellar_Type(2)"][()]

                DELAY_TIMES = fdata['BSE_Double_Compact_Objects']["Coalescence_Time"][()] # Time from ZAMS birth to merger
                Z = np.full(len(M1_CO), met)
                
                fdata.close()

            except Exception:
                print("No DCOs were found")
                
                            
            Z_ALL.append(met)
            
            M1_ZAMS_ALL.append(M1_ZAMS)
            M2_ZAMS_ALL.append(M2_ZAMS)
            
            M1_CO_ALL.append(M1_CO)
            M2_CO_ALL.append(M2_CO)
            
            SMA_CO_ALL.append(SMA_CO)
            ECC_CO_ALL.append(ECC_CO)
            
            TYPE1_CO_ALL.append(TYPE1_CO)
            TYPE2_CO_ALL.append(TYPE2_CO)
            DELAY_TIMES_ALL.append(DELAY_TIMES)
            
            
    return Z_ALL, M1_ZAMS_ALL, M2_ZAMS_ALL, M1_CO_ALL, M2_CO_ALL, \
                SMA_CO_ALL, ECC_CO_ALL, TYPE1_CO_ALL, TYPE2_CO_ALL, DELAY_TIMES_ALL 



def get_sfm(n_systems, metallicities, M1_ZAMS_TOTAL, M2_ZAMS_TOTAL, SFM_CORR=3.671, output_loc=None):
    '''
    Computes the total Star Forming Mass in a galaxy where n_systems were simulated from IMF = Kroupa with [5, 150] M_sol
    SFM_CORR is computed in imf_correction notebook
    '''
    SFM_ALL = np.zeros(len(n_systems))
    for i in range(len(n_systems)):
        index = i*len(metallicities)
        SFM_ALL[i] = (M1_ZAMS_TOTAL[index]+M2_ZAMS_TOTAL[index]) * SFM_CORR
        print(f"Star Forming Mass in {n_systems[i]:.2e} stars: {SFM_ALL[i]:.2e} M_sol \n")
    if output_loc is not None:
        with open(output_loc,'w') as myfile:
            json.dump(list(SFM_ALL), myfile)
      
    return SFM_ALL


def interpolate_n_bbh(sfm, metallicities, M1_BBH_ALL):
    '''
    Returns a 2d interpolator to get the number of BBHs for any given star forming mass and metallicity
    '''
    N_BBH = np.zeros(len(M1_BBH_ALL))
    for i in range(len(N_BBH)):
        N_BBH[i] = len(M1_BBH_ALL[i])
        
    N_BBH_GRID = np.reshape(N_BBH, (len(sfm), len(metallicities)))
    np.shape(N_BBH_GRID)
    
    met_sfm_to_n_bbh = interpolate.interp2d(metallicities, sfm, N_BBH_GRID)
    return met_sfm_to_n_bbh


def met_to_met_weights(met, metallicities):
    '''
    Computes a normal distrbution of metallicities centered around 'met' according to the 
    distribution specified in Fig. 10 and Eq. (48) of https://arxiv.org/pdf/2109.10352.pdf 
    (the COMPAS doc paper).
    '''
    SIGMA_Z = 0.39
    log_metallicities = np.log10(metallicities)
    log_met = np.log10(met)
    weights = np.zeros(len(log_metallicities))
    
    for i in range(len(log_metallicities)):
        mu = log_metallicities[i]
        weights[i] = stats.norm.pdf(log_met, mu, SIGMA_Z)
    weights = weights/np.sum(weights)
    return weights    


def flatten_list_by_met(sfm, metallicities, LIST): 
    '''
    Groups a simulated list of binary parameters (e.g. M1_BBH) by metallicity, flattening across 'n_systems'
    
    Parameters:
        sfm:              list of floats        
            list of star forming masses represented in the simulation
            used as a proxy for n_systems
            
        metallicities:    list of floats
            list of metallicities represented in the simulation
            
        LIST:             list of lists             
            list of Binary properties such that each simulation in a separate list, 
            for a total of len(sfm) x len(metallicities) lists
            
        
    Output: 
        list_complete:     list of lists
            list of Binary properties grouped by only metallicity, i.e. flattened by SFM

    '''
    list_complete = []
    for m in range(len(metallicities)):
        # combine all the binaries for each metallicity
        list_flat = []
        for n in range(len(sfm)):
            l_temp = LIST[get_index(n,m, metallicities)]
            list_flat = [*list_flat, *l_temp]
            
        list_complete.append(list_flat)
    return list_complete


def draw_bbh_from_met(n_met_samples, metallicities, M1_BBH_METS, M2_BBH_METS, DELAY_TIMES_METS):
    '''
    Draws a set of BBHs including their m1, m2, and delay times based on a specified number of samples by metallicity.
    
    Parameters:
        n_met_samples:           list of ints        
            number of BBHs to sample from each metallicity in the 'metallicities' array
            
        metallicities:           list of floats
            list of metallicities that the 'n_met_samples' correspond to
            
        M1_BBH_METS:             list of lists             
            list of all the BBH Primary masses for each metallicity in the 'metallicities' array
            
        M2_BBH_METS:             list of lists             
            list of all the BBH Secondary masses for each metallicity in the 'metallicities' array
            
        DELAY_TIMES_METS:        list of lists             
            list of all the BBH Delay Times for each metallicity in the 'metallicities' array
        
        
    Output: 
       m1_sampled:                list of lists
           sampled primary masses in separate lists corresponding to each metallicity in the 'metallicities' array
       
       m2_sampled:                list of lists
            sampled secondary masses in separate lists corresponding to each metallicity in the 'metallicities' array
       
       delay_time_sampled:        list of lists
           sampled delay times in separate lists corresponding to each metallicity in the 'metallicities' array
    '''
    assert len(metallicities)==len(n_met_samples), "The length of sample array must be the same as the metallicity array"
    
    m1_sampled = []
    m2_sampled = []
    delay_time_sampled = []
    
    for met in range(len(n_met_samples)):
        indices = np.zeros(n_met_samples[met])
        for j in range(len(indices)):
            indices[j] = np.random.randint(0, len(M1_BBH_METS[met]))
        indices = indices.astype(int)
        
        m1_sampled.append(np.array(M1_BBH_METS[met])[indices])
        m2_sampled.append(np.array(M2_BBH_METS[met])[indices])
        delay_time_sampled.append(np.array(DELAY_TIMES_METS[met])[indices])
        
    return m1_sampled, m2_sampled, delay_time_sampled
        
    
def sample_bbh_from_sfm_met(SFM, MET, sfm, metallicities, met_sfm_to_n_bbh, M1_BBH_METS, M2_BBH_METS, DELAY_TIMES_METS):
    '''
    THE MAIN FUNCTION of this file (!)
    Draws a set of BBHs including their m1, m2, and delay times based on a specified metallicity and star forming mass
    
    Parameters:
        SFM:                     float
            The Star Forming Mass for which the BBHs need to be drawn
            
        MET:                     float
            The metallicity for which the BBHs need to be drawn
            
        met_sfm_to_n_bbh:        scipy.interpolate.interpolate.interp2d       
            2d interpolator that can generate the number of BBHs from metallicity and sfm
            this is basically a wrapper for interpolate_n_bbh()
            
        metallicities:           list of floats
            list of metallicities that the 'n_met_samples' correspond to
            
        M1_BBH_METS:             list of lists             
            list of all the BBH Primary masses for each metallicity in the 'metallicities' array
            
        M2_BBH_METS:             list of lists             
            list of all the BBH Secondary masses for each metallicity in the 'metallicities' array
            
        DELAY_TIMES_METS:        list of lists             
            list of all the BBH Delay Times for each metallicity in the 'metallicities' array
        
        
    Output: 
       M1_SAMPLED:                list of floats
           sampled list of primary masses 
       
       M2_SAMPLED:                list of floats
            sampled list of secondary masses
       
       DELAY_TIME_SAMPLED:        list of floats
           sampled list of delay times
    '''
    weights = met_to_met_weights(MET, metallicities)
    n_bbh = met_sfm_to_n_bbh(MET, SFM)
    n_met_samples = [round(num) for num in (weights*n_bbh)]
    
    M1_SAMPLED, M2_SAMPLED, DELAY_TIME_SAMPLED = draw_bbh_from_met(n_met_samples, metallicities, M1_BBH_METS, M2_BBH_METS, DELAY_TIMES_METS)
    
    M1_SAMPLED = flatten(M1_SAMPLED)
    M2_SAMPLED = flatten(M2_SAMPLED)
    DELAY_TIME_SAMPLED = flatten(DELAY_TIME_SAMPLED)
    
    return M1_SAMPLED, M2_SAMPLED, DELAY_TIME_SAMPLED


def flatten(xss):
    return [x for xs in xss for x in xs]
