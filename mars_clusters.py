import fcm
import fcm.atmosphere as atm
import matplotlib.pyplot as plt
import numpy as np
from random import randint
import pandas as pd
import time
import argparse
import gtools as gt
from gtools import crater_tools

# Optional arguments
parser = argparse.ArgumentParser(description='Program for simulating distribution of small clusters/craters on Mars')
parser.add_argument('--id', '-I', action="store", dest="model_id", required=True)
parser.add_argument('--restart', '-R', action="store_true", dest="restart", default=False)
parser.add_argument('--verbose', '-V', action="store_true", dest="verbose", default=False)
parser.add_argument('--craters', '-C', action="store_true", dest="save_craters", default=False)
parser.add_argument('--scaling', '-S', action="store", dest="cscaling", default="CS")
parser.add_argument('--lift', '-L', action="store", dest="lift_coef", default=0.02, type=float)
parser.add_argument('--fragments', '-F', action="store_true", dest="save_frag", default=False)
parser.add_argument('--velocity', '-v', action="store", dest="velocity", default=-1., type=float)
parser.add_argument('--angle', '-t', action="store", dest="angle", default=-1., type=float)
parser.add_argument('--atmos', '-A', action="store", dest="atmos", default=1, type=int)
parser.add_argument('--break', '-B', action="store", dest="breakm", default=2, type=int)
parser.add_argument('--alpha', '-a', action="store", dest="alpha", default=0.926, type=float)
parser.add_argument('--minmass', '-m', action="store", dest="minmass", default=15., type=float)
parser.add_argument('--dlo', '-l', action="store", dest="dlo", default=1400., type=float)
parser.add_argument('--dup', '-u', action="store", dest="dup", default=4000., type=float)
parser.add_argument('--med_strength', '-s', action="store", dest="med_strength", default=330., type=float)
parser.add_argument('--width', '-w', action="store", dest="width", default=1.25, type=float)
parser.add_argument('--minrad', '-c', action="store", dest="crat_min", default=0.5, type=float)
parser.add_argument('--ablation', '-e', action="store", dest="cab", default=4.2E-8, type=float)
parser.add_argument('--drag', '-d', action="store", dest="cd", default=1., type=float)
parser.add_argument('--ss_disp', '-x', action="store", dest="ss_disp", default=0.5, type=float)
parser.add_argument('--fm_disp', '-y', action="store", dest="fm_disp", default=0.9, type=float)
parser.add_argument('--fvcmin', '-f', action="store", dest="fvc_min", default=1.1, type=float)
parser.add_argument('--fvcmax', '-g', action="store", dest="fvc_max", default=1.1, type=float)
parser.add_argument('--semin', '-j', action="store", dest="se_min", default=0.25, type=float)
parser.add_argument('--semax', '-k', action="store", dest="se_max", default=0.25, type=float)
parser.add_argument('--impacts', '-i', action="store", dest="impacts", default=10000, type=int)
parser.add_argument('--numc', '-n', action="store", dest="numc", default=200, type=int)

def compile_param_dict(options):
    
    hp = {}
    hp['Samples'] = n
    hp['Atmosphere'] = options.atmos
    hp['Mass SFD exp.'] = options.alpha
    hp['Min. mass'] = options.minmass
    hp['Min. density'] = options.dlo
    hp['Max. density'] = options.dup
    hp['Med. strength'] = options.med_strength
    hp['Var. strength'] = options.width
    hp['Strength exp. (min)'] = options.se_min
    hp['Strength exp. (max)'] = options.se_max
    hp['Frag sep coef (min)'] = options.fvc_min
    hp['Frag sep coef (max)'] = options.fvc_max
    hp['Crater rad. (min)'] = options.crat_min
    hp['Ablation coef.'] = options.cab
    hp['Drag coef.'] = options.cd
    hp['Crater scaling'] = options.cscaling
    hp['Strength scaling disp.'] = options.ss_disp
    hp['Fragment mass disp.'] = options.fm_disp
    hp['Break-up mode'] = options.breakm
    hp['Lift coef.'] = options.lift_coef

    return hp  

# Model options
options = parser.parse_args()
model_id = options.model_id          # Model ID for record keeping
n=100000                             # Samples in distributions to draw from
atmos = options.atmos                # With or without atmosphere
alpha = options.alpha                # Power-law exponent in mass distribution
minmass = options.minmass            # Minimum mass
dup = options.dup                    # Density distribution bounds
dlo = options.dlo    
med_strength = options.med_strength  # Median strength
width = options.width                # Strength variation
if width >= 0.:
    lognormalstrength=False          # width is half-width of log-uniform dist.
else:
    lognormalstrength=True           # width is 3-sigma of log-normal dist.
    width = -width
fvc_min = options.fvc_min            # Frag. sep. coef. (min)
fvc_max = options.fvc_max            # Frag. sep. coef. (max)
se_min = options.se_min              # Strength exponent (min)
se_max = options.se_max              # Strength exponent (max)
crater_radius_min = options.crat_min # Minimum crater radius to track in model
ablation_parameter = options.cab     # Ablation coefficient
drag_coef = options.cd
lift_coef = options.lift_coef
fm_disp = options.fm_disp
ss_disp = options.ss_disp
break_mode = options.breakm
impacts = options.impacts            # Max samples in sim.
num_break = options.numc             # Target number of D_eff > 10 craters
if options.cscaling == 'HS':
    c_scaling = 'hard_soil'
elif options.cscaling == 'DS':
    c_scaling = 'cohesionless_material'
elif options.cscaling == 'YS':
    c_scaling = 'dry_soil'
else:
    c_scaling = 'cohesive_soil'


# Save hyper parameters in a dictionary (needs fixing for FCM)
hyper_params = compile_param_dict(options)

# Record the hyper paramters. . .
record = pd.DataFrame()
#record = pd.read_csv('record.csv', index_col='ID')
sample = pd.Series(hyper_params, name=model_id)
record = record.append(sample)
record.to_csv('record.csv')

# Define output file names
outdir   = './'  #Specify output directory here
outfile  = outdir+'output-'+model_id+'.csv' 
cratfile = outdir+'craters-'+model_id+'.csv'
fragfile = outdir+'fragments-'+model_id+'.csv'

# Load atmospheric density vs elevation data
if atmos:
    atmosphere = atm.static_martian_atmosphere()
else:
    atmosphere = atm.exponential(rho0=1E-9, hmax=1E6, scale_height=1E6)


# Define the parameter distributions

# Velocity
if options.velocity < 0.:
    container = np.load('gtools/velocity_distribution.npz')
    values = container["velocity"]
    cdf = container["cdf"]
    container.close()
    seed = 23456
    generator = np.random.default_rng(seed)
    uniform_samples = generator.uniform(size=n)
    indices = np.searchsorted(cdf, uniform_samples)
    velocity = values[indices]
else:
    velocity = [options.velocity]*n
    
# Angle
if options.angle < 0.:
    P = np.random.random(n)
    angle = np.rad2deg(np.arcsin(np.sqrt(P)))
else:
    angle = [options.angle]*n

# Density, ablation
density = np.random.uniform(dlo, dup, n)
ablation = np.random.uniform(1.E-8, ablation_parameter, n)

# strength
if lognormalstrength:
    strength = med_strength * 10**np.random.normal(scale=width/3, size=n)
else:
    strength = 10**np.random.uniform(np.log10(med_strength)-width,
                                        np.log10(med_strength)+width,n)
      
# mass - used together with density to give radius
if alpha > 0:
    generator = np.random.default_rng(123)
    mass = minmass * (generator.pareto(alpha, n) + 1)
else:
    mass = [minmass]*n
    
# strength exponent
ses = np.random.uniform(se_min,se_max,n)

# separation coefficient
fvcs = np.random.uniform(fvc_min,fvc_max,n)

# Simulate impacts
if options.restart:
    otmp = pd.read_csv(outfile)
    istart = len(otmp)
    otmp = otmp[otmp['Effective Diameter [m]'] > 10]
    num_g10 = len(otmp)
    num_sg10 = len(otmp[otmp['No. of Craters'] == 1])
    nosave = False
    del otmp
    print("Restarting run: ",model_id," from sample: ",istart)
    print("Number of large craters so far: ", num_g10, num_sg10)
else:
    istart = 0
    num_g10 = 0              # Number of craters with D_eff > 10
    num_sg10 = 0             # Number of ind. craters with D_eff > 10
    nosave = True
output = pd.DataFrame()  # Output from one sample
inputs = {}              # Inputs from one sample
ii = 0; jj = 0           # For counting number of samples since last save
for i in range(istart,impacts):
    
    # Random selection of parameters
    inputs['Velocity [km/s]']=np.random.choice(velocity)
    inputs['Angle']=np.random.choice(angle)
    
    # Density strength and ablation coefficient from single dist., or dependent on type
    d = np.random.choice(density)
    Y = np.random.choice(strength)
    C_ab = np.random.choice(ablation)
    inputs['Ablation coef.']=C_ab
    inputs['Density [kg/m3]']=d 

    # Impactor density, mass and radius
    m=np.random.choice(mass)
    if (m > 20E3 * minmass):  # Limit the largest mass that is considered for efficiency
        m=np.random.choice(mass)
    inputs['Radius [m]']=np.cbrt(3*m/(4*np.pi*d))
    inputs['Mass [kg]']=m
  
    # Other coefficients
    inputs['Strength exponent']=np.random.choice(ses)
    inputs['Frag. sep. coef.']=np.random.choice(fvcs)
    
    # Model parameters
    parameters = fcm.FCMparameters(g0=3.72, Rp=3390, atmospheric_density=atmosphere, 
                                   frag_velocity_coeff=inputs['Frag. sep. coef.'],
                                   ablation_coeff=C_ab, drag_coeff=drag_coef, lift_coeff=lift_coef,
                                   min_crater_radius=crater_radius_min,
                                   strengh_scaling_disp=ss_disp, fragment_mass_disp=fm_disp, 
                                   cratering_params=c_scaling)
           
    # Define a meteoroid
    inputs['Strength [kPa]']=max(0.1, min(Y, parameters.max_strength))
    impactor = fcm.FragmentationMeteoroid(velocity=inputs['Velocity [km/s]'], angle=inputs['Angle'], 
                                              density=inputs['Density [kg/m3]'], radius=inputs['Radius [m]'], 
                                              strength=inputs['Strength [kPa]'], 
                                              strength_scaler=inputs['Strength exponent'],
                                              fragment_mass_fractions=gt.random_mass_fractions(break_mode))
        
    # Run the fragmentation model
    start = time.time()
    if options.save_frag:
        results = fcm.simulate_impact(parameters, impactor, h_start=100, seed=randint(0,100), final_states=True)
    else:
        results = fcm.simulate_impact(parameters, impactor, h_start=100, seed=randint(0,100))
    end = time.time()
    inputs['CPU Time (s)'] = end - start

    # Process and record the data
    characteristics = crater_tools.cluster_characteristics(results.craters)
    sample = pd.Series({**inputs, **characteristics}, name=str(i))
    output = output.append(sample)

    # If desired, save the crater information and/or final states to compendium
    if options.save_frag and results.final_states is not None:
        fs = results.final_states.reset_index(drop=True)
        tdf = fs.assign(ID=str(i)).set_index('ID',append=True).swaplevel(0,1)
        if jj == 0:
            fss = tdf
        else:
            fss = fss.append(tdf)
        jj += 1
            
    if options.save_craters and results.craters is not None:
        try:
            cs = results.craters.drop(columns=['IDs'])
        except:
            cs = results.craters
        tdf = cs.assign(ID=str(i)).set_index('ID',append=True).swaplevel(0,1)
        if ii == 0:
            css = tdf
        else:
            css = css.append(tdf)
        ii += 1
    
    # Monitor number of craters / clusters with D_eff > 10 as stopping criterion
    if characteristics['Effective Diameter [m]'] > 10:
        num_g10 += 1
        if characteristics['No. of Craters'] == 1:
            num_sg10 += 1

    # Periodically save progress and print status
    if (i+1) % 20 == 0:

        # Saving data
        if nosave:
            output.to_csv(outfile)
            if options.save_frag:
                fss.to_csv(fragfile)
                jj = 0
            if options.save_craters:
                css.to_csv(cratfile)
                ii = 0
            output = pd.DataFrame()
            nosave = False
        else:
            output.to_csv(outfile, mode='a', header=False)
            output = pd.DataFrame()
            if options.save_frag:
                fss.to_csv(fragfile, mode='a', header=False)
                jj = 0
            if options.save_craters:
                css.to_csv(cratfile, mode='a', header=False)
                ii = 0
            
        # Print status in verbose mode
        if options.verbose:
            try:
                print(i+1, num_g10, num_sg10, num_g10/num_sg10-1.)
            except:
                print(i+1, num_g10, num_sg10)

    # Early stopping criterion: desired number of large clusters
    if num_g10 == num_break:
        break
   
# Save final datasets
if nosave:
    output.to_csv(outfile)
    if options.save_frag:
        try:
            fss.to_csv(fragfile)
        except:
            print('Could not save fragment file')
    if options.save_craters:
        try:
            css.to_csv(cratfile)
        except:
            print('Could not save craters file')
else:
    output.to_csv(outfile, mode='a', header=False)
    if options.save_frag:
        fss.to_csv(fragfile, mode='a', header=False)
    if options.save_craters:
        css.to_csv(cratfile, mode='a', header=False)
