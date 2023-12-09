from PlanetClass import Planet
from SolarSystemClass import SolarSystemPP
import numpy as np
from planetary_data import earth, sun, jupiter
import matplotlib.pyplot as plt
import numba as nb


MSun = 1.988e30
ME = 5.972e24
AU = 1.496e8 # km
G = 6.67430e-20  # km**3/ kg / s**2

Sun = Planet('Sun', sun['radius'], 1.3271244004193938E+11, 'red', 1., 1., MSun, 0., 0., 0., 0.)
Earth = Planet('Earth', earth['radius'], 5.972e24 * G, 'blue',1., 1., ME, AU, 0., 0., 30.)
Jupiter = Planet('Jupiter', jupiter['radius'], 1.26686e8, 'red', 1., 1., 1.898e27, 5.20336301 * AU, 0., 0., 13.06)
Mercury = Planet("Mercury", 2439.7, 0.3301e24 * G, "gray", 57.909e6, 0.02056,  0.3301e24, -2.21207e7, -6.68243e7, 3.66622e1, -1.23026e1)

def rr(x, y):
    return np.sqrt(x**2 + y**2)



# Euler_results
nsteps = np.logspace(3.5, 6.5, 10, dtype= int)
apx = np.array([])
apy = np.array([])
perix = np.array([])
periy = np.array([])

for steps in nsteps:
    # initiate space
    SS = SolarSystemPP(1)
    SS.addBody(Sun)
    SS.addBody(Mercury)

    # start orbiting, 1 years
    SS.Orbit(steps, 3.1536e7, 'Euler', GR = True)

    # Store mercury data
    M_x = Mercury.trajectory[:,0]
    M_y = Mercury.trajectory[:,1]
    M_r = rr(M_x, M_y)

    # pick out only the peri/aphelion point to store
    M_peri_x = np.array([])
    M_peri_y = np.array([])
    M_ap_x = np.array([])
    M_ap_y = np.array([])
    for i in range(len(M_x) - 2):
        if M_r[i+1] > M_r[i+2] and M_r[i] < M_r[i+1]:
            M_peri_x = np.insert(M_peri_x, 0, M_x[i+1], axis=0)
            M_peri_y = np.insert(M_peri_y, 0, M_y[i+1], axis=0)

        elif M_r[i+1] < M_r[i+2] and M_r[i] > M_r[i+1]:
            M_ap_x = np.insert(M_ap_x, 0, M_x[i+1], axis=0)
            M_ap_y = np.insert(M_ap_y, 0, M_y[i+1], axis=0)


    # further pick out the initial/final peri/aphelion
    apx = np.insert(apx, 0, M_ap_x[-1], axis=0)
    apx = np.insert(apx, 0, M_ap_x[0], axis=0)

    apy = np.insert(apy, 0, M_ap_y[-1], axis=0)
    apy = np.insert(apy, 0, M_ap_y[0], axis=0)

    perix = np.insert(perix, 0, M_peri_x[-1], axis=0)
    perix = np.insert(perix, 0, M_peri_x[0], axis=0)

    periy = np.insert(periy, 0, M_peri_y[-1], axis=0)
    periy = np.insert(periy, 0, M_peri_y[0], axis=0)

# reverse the array to get right order
apx = apx[::-1]
apy = apy[::-1]
perix = perix[::-1]
periy = periy[::-1]
nsteps = np.pad(nsteps, (0,10), "constant", constant_values=(0, 0))
# save the end results of initial/final peri/aphelion data
# for further analysis 
np.savetxt('M_Euler_results.csv', (nsteps, apx, apy, perix, periy), delimiter=',')




# RK4_results
tol_list = np.logspace(-1, 1.7, 10)

apx = np.array([])
apy = np.array([])
perix = np.array([])
periy = np.array([])

numsteps = np.array([])

for toll in tol_list:
    # initiate space
    SS = SolarSystemPP(1)
    SS.addBody(Sun)
    SS.addBody(Mercury)

    # start orbiting, 2 years
    SS.Orbit_RK_GR(3.1536e7, tol = toll)

    # Store mercury data
    M_x = Mercury.trajectory[:,0]
    M_y = Mercury.trajectory[:,1]
    M_r = rr(M_x, M_y)

    # pick out only the peri/aphelion point to store
    M_peri_x = np.array([])
    M_peri_y = np.array([])
    M_ap_x = np.array([])
    M_ap_y = np.array([])
    for i in range(len(M_x) - 2):
        if M_r[i+1] > M_r[i+2] and M_r[i] < M_r[i+1]:
            M_peri_x = np.insert(M_peri_x, 0, M_x[i+1], axis=0)
            M_peri_y = np.insert(M_peri_y, 0, M_y[i+1], axis=0)

        elif M_r[i+1] < M_r[i+2] and M_r[i] > M_r[i+1]:
            M_ap_x = np.insert(M_ap_x, 0, M_x[i+1], axis=0)
            M_ap_y = np.insert(M_ap_y, 0, M_y[i+1], axis=0)


    # further pick out the initial/final peri/aphelion
    apx = np.insert(apx, 0, M_ap_x[-1], axis=0)
    apx = np.insert(apx, 0, M_ap_x[0], axis=0)

    apy = np.insert(apy, 0, M_ap_y[-1], axis=0)
    apy = np.insert(apy, 0, M_ap_y[0], axis=0)

    perix = np.insert(perix, 0, M_peri_x[-1], axis=0)
    perix = np.insert(perix, 0, M_peri_x[0], axis=0)

    periy = np.insert(periy, 0, M_peri_y[-1], axis=0)
    periy = np.insert(periy, 0, M_peri_y[0], axis=0)

    numsteps = np.insert(numsteps, 0, len(M_x), axis = 0)

# reverse the array to get right order
apx = apx[::-1]
apy = apy[::-1]
perix = perix[::-1]
periy = periy[::-1]
numsteps = numsteps[::-1]
numsteps = np.pad(numsteps, (0,10), "constant", constant_values=(0, 0))
# save the end results of initial/final peri/aphelion data
# for further analysis 
np.savetxt('M_RK4_results.csv', (numsteps, apx, apy, perix, periy), delimiter=',')


