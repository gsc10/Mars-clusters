import numpy as np

# Functions
def choose_impactor_w_comets(x):
    
    if x < 0.05:
        d_i = 7500.*10**np.random.normal(0.0, 0.01)
        C_ab_i = 7.E-8
        Y_i = 1.E4*10**np.random.normal(0.0, 0.2)
    elif x < .35:
        d_i = 3400.*np.random.normal(1.0, 0.15)
        C_ab_i = 1.4E-8
        Y_i = 1.E3*10**np.random.normal(0., 0.4)
    elif x < 0.70:
        d_i = 2600.*np.random.normal(1.0, 0.15)
        C_ab_i = 4.2E-8
        Y_i = 3.E2*10**np.random.normal(0.0, 0.4)
    else:
        d_i = 1000.*10**np.random.normal(0.0, 0.2)
        C_ab_i = 1.E-7 #2.1E-7
        Y_i = 3.E1*10**np.random.normal(0.0, 0.2)
    
    return d_i, C_ab_i, Y_i

def choose_impactor(x, Y_oc_med=1.E3):
    
    if x < 0.05:
        # Irons
        d_i = 7500.*10**np.random.normal(0.0, 0.01)
        C_ab_i = 7.E-8
        Y_i = 1.E4*10**np.random.normal(0.0, 0.2)
    elif x < .45:
        # Ordinary chondrites
        d_i = 3400.*np.random.normal(1.0, 0.15)
        C_ab_i = 1.4E-8
        Y_i = Y_oc_med*10**np.random.normal(0., 0.45)
    else:
        # Carbonaceous chondrites
        d_i = 2600.*np.random.normal(1.0, 0.15)
        C_ab_i = 4.2E-8
        Y_i = 3.E2*10**np.random.normal(0.0, 0.4)
    
    return d_i, C_ab_i, Y_i

def random_mass_fractions(break_mode):
    if break_mode == 1:
        l=np.random.randint(1,9)/10
    else:
        l=1/break_mode
    return (l,1-l)
