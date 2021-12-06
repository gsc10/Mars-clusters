import numpy as np

def strength_scaling_momentum(Dr, rho=1500., delta=2000., vvert=1000., strength=1.E5):
    return 821.*(Dr)**3 * (strength/1E5)**0.62 * (rho/1500.)**0.59 \
                        * (delta/2000.)**-0.2 * (vvert/1000.)**-0.23

def gravity_scaling_momentum(Dr, rho=1500., delta=2000., vvert=1000., gravity=3.71):
    return 79.4*(Dr)**3.6 * (gravity/3.71)**0.61 * (rho/1500.)**1.19 \
                          * (delta/2000.)**-0.19 * (vvert/1000.)**-0.23

def holsapple(radius, grav, rho_t, rho_i, U, Y0, K1, K2, mu, nu, Kr, Kd):
    # impactor mass (assumes a sphere)
    m = 4. / 3. * np.pi * rho_i * radius**3
    # gravity-scaled impactor size (\pi_2)
    pi2 = radius * grav / U**2
    # strength-scaled impactor size (\pi_3)
    pi3 = Y0 / (rho_t * U**2)
    # target-impactor density ratio
    pi4 = rho_t / rho_i
    # exponents
    exp1 = (6 * nu - 2. - mu) / (3. * mu)
    exp2 = (6 * nu - 2.) / (3. * mu)
    exp3 = (2 + mu) / 2.
    exp4 = -3. * mu / (2 + mu)
    # Cratering efficiency (Volume)
    piV = K1 * (pi2 * pi4**exp1 + K2 * (pi3 * pi4**exp2)**exp3)**exp4
    # Crater volume (below preimpact surface)
    V = piV * m / rho_t
    # Crater radius (at preimpact level)
    Rt = Kr * V**(1. / 3.)
    # Crater rim radius
    Rf = 1.3 * Rt
    # scaled crater diameter (rim)
    piD = 2.*Rf*(rho_t/m)**(1./3.)
    # Crater depth (floor to preimpact level)
    d = Kd / Kr * Rt
    # Kinetic energy of impactor
    E = 0.5 * m * U**2
    # Momentum of impactor
    L = m * U
    return pi2, pi3, piV, piD, V, Rt, Rf, d, E, L
