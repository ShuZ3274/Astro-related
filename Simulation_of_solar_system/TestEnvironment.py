# -*- coding: utf-8 -*-
# Custom classes and methods
from PlanetClass import Body, Planet
from SolarSystemClass import SolarSystem, SolarSystemPP
from BHClass import SolarSystemBH, displaySystem

# Data for known/identified celstial objects
from planetary_data import earth, sun, jupiter, mars, venus, mercury
from CelestialData import Earth, Jupiter, Mars, Mercury, Venus

# Standard libraries
import matplotlib.pyplot as plt
import numpy as np
import math
import time

from mpl_toolkits.mplot3d import Axes3D

MSun = 1.988500e30
ME = 5.972e24
AU = 1.496e8 # km
G_meters = 6.67430e-11       # m**3 / kg / s**2
G        = G_meters * 10**-9 # km**3/ kg / s**2

'''
Sun =       Planet('Sun',       sun['radius'],      1.3271244004193938E+11, 'red',      1.,          1.,        MSun,               0.,             0.,             0.,             0.)
Earth =     Planet('Earth',     earth['radius'],    5.972e24 * G,           'blue',     1.,          1.,        earth['mass'],      earth['x'],     earth['y'],     earth['u'],     earth['v'])
Jupiter =   Planet('Jupiter',   jupiter['radius'],  1.26686e8,              'red',      1.,          1.,        jupiter['mass'],    jupiter['x'],   jupiter['y'],   jupiter['u'],   jupiter['v'])
Mercury =   Planet("Mercury",   2439.7,             0.3301e24 * G,          "gray",     57.909e6,    0.02056,   mercury['mass'],    mercury['x'],   mercury['y'],   mercury['u'],   mercury['v'])
Mars =      Planet('Mars',      mars['radius'],     mars['mu'],             'green',    1.,          1.,        mars['mass'],       mars['x'],      mars['y'],      mars['u'],      mars['v'])
Venus =     Planet('Venus',     venus['radius'],    venus['mu'],            'green',    1.,          1.,        venus['mass'],      venus['x'],     venus['y'],     venus['u'],     venus['v'])


SS = SolarSystemBH(15 * AU, 0.5)


SS.addBody(Sun)
SS.addBody(Earth)
SS.addBody(Mercury)
SS.addBody(Jupiter)
SS.addBody(Mars)
SS.addBody(Venus)

SS.Orbit(12*365, 12 * 365 * 24 * 3600)

plt.figure()
fig, ax = plt.subplots()
ax.plot(SS.getBodies()[0].trajectory[:, 0], SS.getBodies()[0].trajectory[:, 1], color='yellow', ms = 50, label='Sun')
ax.plot(SS.getBodies()[2].trajectory[:, 0], SS.getBodies()[2].trajectory[:, 1], color='brown', ms = 10, label='Mercury')
ax.plot(SS.getBodies()[5].trajectory[:, 0], SS.getBodies()[5].trajectory[:, 1], color='gray', ms = 20, label='Venus')
ax.plot(SS.getBodies()[1].trajectory[:, 0], SS.getBodies()[1].trajectory[:, 1], color='green', ms = 20, label='Earth')
ax.plot(SS.getBodies()[4].trajectory[:, 0], SS.getBodies()[4].trajectory[:, 1], color='red', ms = 15, label='Mars')
ax.plot(SS.getBodies()[3].trajectory[:, 0], SS.getBodies()[3].trajectory[:, 1], color='purple', ms = 35, label='Jupiter')
plt.legend(loc = 'upper right')
ax.set_xlabel('x[km]')
ax.set_ylabel('y[km]')
plt.show()


for body in SS.getBodies():
    plt.plot(body.trajectory[:, 0], body.trajectory[:, 1])
plt.show()


pos1 = []
    '''

SS = SolarSystemBH(15*AU, 0)    

rng1 = np.random.default_rng(seed = 82636473836)
rng2 = np.random.default_rng(seed = 84634357238)

Sun =       Planet('Sun',       sun['radius'],      1.3271244004193938E+11, 'red',      1.,          1.,        MSun,               0.,             0.,             0.,             0.)
Earth =     Planet('Earth',     earth['radius'],    5.972e24 * G,           'blue',     1.,          1.,        earth['mass'],      earth['x'],     earth['y'],     earth['u'],     earth['v'])
Jupiter =   Planet('Jupiter',   jupiter['radius'],  1.26686e8,              'red',      1.,          1.,        jupiter['mass'],    jupiter['x'],   jupiter['y'],   jupiter['u'],   jupiter['v'])
Mercury =   Planet("Mercury",   2439.7,             0.3301e24 * G,          "gray",     57.909e6,    0.02056,   mercury['mass'],    mercury['x'],   mercury['y'],   mercury['u'],   mercury['v'])
Mars =      Planet('Mars',      mars['radius'],     mars['mu'],             'green',    1.,          1.,        mars['mass'],       mars['x'],      mars['y'],      mars['u'],      mars['v'])
Venus =     Planet('Venus',     venus['radius'],    venus['mu'],            'green',    1.,          1.,        venus['mass'],      venus['x'],     venus['y'],     venus['u'],     venus['v'])

nbodies = 400
rand_r_trojans = (rng1.normal(loc = 0., scale = 0.35, size = nbodies) + 5 ) * AU
rand_theta_trojans = rng1.normal(loc = 0., scale = 0.5, size = nbodies) * np.pi / 6 + 0.63306 # math.atan2(Jupiter['y'], Jupiter['x'])     
rand_theta_trojans[:nbodies // 2] += np.pi / 3
rand_theta_trojans[nbodies // 2:9 * (nbodies // 10)] -= np.pi / 3
rand_theta_trojans[9 * (nbodies // 10):] += np.pi

rand_r_belt = (rng2.normal(loc = 0., scale = 0.5, size = nbodies ) + 2.66 ) * AU
rand_theta_belt = rng2.random(nbodies) * 2 *  np.pi 

x_rand_trojans = rand_r_trojans * np.cos(rand_theta_trojans)
y_rand_trojans = rand_r_trojans * np.sin(rand_theta_trojans)

x_rand_belt = rand_r_belt * np.cos(rand_theta_belt)
y_rand_belt = rand_r_belt * np.sin(rand_theta_belt)

speed_rand_trojans = np.sqrt(G * MSun / rand_r_trojans) 

u_rand = -np.sin(rand_theta_trojans) * speed_rand_trojans
v_rand = np.cos(rand_theta_trojans) * speed_rand_trojans

speed_rand_belt = np.sqrt(G * MSun / rand_r_belt) 

u_rand_belt = -np.sin(rand_theta_belt) * speed_rand_belt
v_rand_belt = np.cos(rand_theta_belt) * speed_rand_belt

mass_rand = 10**(np.random.rand(nbodies) * 18)
 

SS.addBody(Sun)
SS.addBody(Mercury)
SS.addBody(Venus)
SS.addBody(Earth)
SS.addBody(Mars)
SS.addBody(Jupiter)
    

for i in range(nbodies):
    belt = Planet(f'Belt {i}', 1., 1., 'str', 1., 1., 0.1, x_rand_belt[i], y_rand_belt[i], u_rand[i], v_rand[i])
    trojan = Planet(f'Trojan {i}', 1., 1., 'str', 1., 1., 0.1, x_rand_trojans[i], y_rand_trojans[i], u_rand[i], v_rand[i])
    SS.addBody(trojan)
    SS.addBody(belt)
    
SS.buildTree()
displaySystem(SS)
    
''' 
    for i in range(nbodies):
        belt = Planet(f'Belt {i}', 1., 1., 'str', 1., 1., 0.1, x_rand_belt[i], y_rand_belt[i], u_rand[i], v_rand[i])
        trojan = Planet(f'Trojan {i}', 1., 1., 'str', 1., 1., 0.1, x_rand_trojans[i], y_rand_trojans[i], u_rand[i], v_rand[i])
        SS.addBody(trojan)
        SS.addBody(belt)
        
        
    SS.Evolve(3600)
    pos_true = SS.getPositions()
    pos1.append(pos_true)
    
    
pos2 = []
    
for n in np.linspace(1, 3, 7, endpoint = True):
    
    rng1 = np.random.default_rng(seed = 82636473836)
    rng2 = np.random.default_rng(seed = 84634357238)
    
    Sun =       Planet('Sun',       sun['radius'],      1.3271244004193938E+11, 'red',      1.,          1.,        MSun,               0.,             0.,             0.,             0.)
    Earth =     Planet('Earth',     earth['radius'],    5.972e24 * G,           'blue',     1.,          1.,        earth['mass'],      earth['x'],     earth['y'],     earth['u'],     earth['v'])
    Jupiter =   Planet('Jupiter',   jupiter['radius'],  1.26686e8,              'red',      1.,          1.,        jupiter['mass'],    jupiter['x'],   jupiter['y'],   jupiter['u'],   jupiter['v'])
    Mercury =   Planet("Mercury",   2439.7,             0.3301e24 * G,          "gray",     57.909e6,    0.02056,   mercury['mass'],    mercury['x'],   mercury['y'],   mercury['u'],   mercury['v'])
    Mars =      Planet('Mars',      mars['radius'],     mars['mu'],             'green',    1.,          1.,        mars['mass'],       mars['x'],      mars['y'],      mars['u'],      mars['v'])
    Venus =     Planet('Venus',     venus['radius'],    venus['mu'],            'green',    1.,          1.,        venus['mass'],      venus['x'],     venus['y'],     venus['u'],     venus['v'])
    
    nbodies = int(10**n)
    rand_r_trojans = (rng1.normal(loc = 0., scale = 0.35, size = nbodies) + 5 ) * AU
    rand_theta_trojans = rng1.normal(loc = 0., scale = 0.5, size = nbodies) * np.pi / 6 + 0.63306 # math.atan2(Jupiter['y'], Jupiter['x'])     
    rand_theta_trojans[:nbodies // 2] += np.pi / 3
    rand_theta_trojans[nbodies // 2:9 * (nbodies // 10)] -= np.pi / 3
    rand_theta_trojans[9 * (nbodies // 10):] += np.pi
    
    rand_r_belt = (rng2.normal(loc = 0., scale = 0.5, size = nbodies ) + 2.66 ) * AU
    rand_theta_belt = rng2.random(nbodies) * 2 *  np.pi 
    
    x_rand_trojans = rand_r_trojans * np.cos(rand_theta_trojans)
    y_rand_trojans = rand_r_trojans * np.sin(rand_theta_trojans)
    
    x_rand_belt = rand_r_belt * np.cos(rand_theta_belt)
    y_rand_belt = rand_r_belt * np.sin(rand_theta_belt)
    
    speed_rand_trojans = np.sqrt(G * MSun / rand_r_trojans) 
    
    u_rand = -np.sin(rand_theta_trojans) * speed_rand_trojans
    v_rand = np.cos(rand_theta_trojans) * speed_rand_trojans
    
    speed_rand_belt = np.sqrt(G * MSun / rand_r_belt) 
    
    u_rand_belt = -np.sin(rand_theta_belt) * speed_rand_belt
    v_rand_belt = np.cos(rand_theta_belt) * speed_rand_belt
    
    mass_rand = 10**(np.random.rand(nbodies) * 18)
    
    SS = SolarSystemPMBH2(15 * AU, 0)   
    
    SS.addBody(Sun)
    SS.addBody(Earth)
    SS.addBody(Mercury)
    SS.addBody(Jupiter)
    SS.addBody(Mars)
    SS.addBody(Venus)
    
    
    
    for i in range(nbodies):
        belt = Planet(f'Belt {i}', 1., 1., 'str', 1., 1., 0.1, x_rand_belt[i], y_rand_belt[i], u_rand[i], v_rand[i])
        trojan = Planet(f'Trojan {i}', 1., 1., 'str', 1., 1., 0.1, x_rand_trojans[i], y_rand_trojans[i], u_rand[i], v_rand[i])
        SS.addBody(trojan)
        SS.addBody(belt)
        
        
    SS.Evolve(3600, method = 'Euler')
    pos = SS.getPositions()
    pos2.append(pos)
    
pos3 = []
    
for n in np.linspace(1, 3, 7, endpoint = True):
    
    rng1 = np.random.default_rng(seed = 82636473836)
    rng2 = np.random.default_rng(seed = 84634357238)
    
    Sun =       Planet('Sun',       sun['radius'],      1.3271244004193938E+11, 'red',      1.,          1.,        MSun,               0.,             0.,             0.,             0.)
    Earth =     Planet('Earth',     earth['radius'],    5.972e24 * G,           'blue',     1.,          1.,        earth['mass'],      earth['x'],     earth['y'],     earth['u'],     earth['v'])
    Jupiter =   Planet('Jupiter',   jupiter['radius'],  1.26686e8,              'red',      1.,          1.,        jupiter['mass'],    jupiter['x'],   jupiter['y'],   jupiter['u'],   jupiter['v'])
    Mercury =   Planet("Mercury",   2439.7,             0.3301e24 * G,          "gray",     57.909e6,    0.02056,   mercury['mass'],    mercury['x'],   mercury['y'],   mercury['u'],   mercury['v'])
    Mars =      Planet('Mars',      mars['radius'],     mars['mu'],             'green',    1.,          1.,        mars['mass'],       mars['x'],      mars['y'],      mars['u'],      mars['v'])
    Venus =     Planet('Venus',     venus['radius'],    venus['mu'],            'green',    1.,          1.,        venus['mass'],      venus['x'],     venus['y'],     venus['u'],     venus['v'])
    
    nbodies = int(10**n)
    rand_r_trojans = (rng1.normal(loc = 0., scale = 0.35, size = nbodies) + 5 ) * AU
    rand_theta_trojans = rng1.normal(loc = 0., scale = 0.5, size = nbodies) * np.pi / 6 + 0.63306 # math.atan2(Jupiter['y'], Jupiter['x'])     
    rand_theta_trojans[:nbodies // 2] += np.pi / 3
    rand_theta_trojans[nbodies // 2:9 * (nbodies // 10)] -= np.pi / 3
    rand_theta_trojans[9 * (nbodies // 10):] += np.pi
    
    rand_r_belt = (rng2.normal(loc = 0., scale = 0.5, size = nbodies ) + 2.66 ) * AU
    rand_theta_belt = rng2.random(nbodies) * 2 *  np.pi 
    
    x_rand_trojans = rand_r_trojans * np.cos(rand_theta_trojans)
    y_rand_trojans = rand_r_trojans * np.sin(rand_theta_trojans)
    
    x_rand_belt = rand_r_belt * np.cos(rand_theta_belt)
    y_rand_belt = rand_r_belt * np.sin(rand_theta_belt)
    
    speed_rand_trojans = np.sqrt(G * MSun / rand_r_trojans) 
    
    u_rand = -np.sin(rand_theta_trojans) * speed_rand_trojans
    v_rand = np.cos(rand_theta_trojans) * speed_rand_trojans
    
    speed_rand_belt = np.sqrt(G * MSun / rand_r_belt) 
    
    u_rand_belt = -np.sin(rand_theta_belt) * speed_rand_belt
    v_rand_belt = np.cos(rand_theta_belt) * speed_rand_belt
    
    mass_rand = 10**(np.random.rand(nbodies) * 18)
    
    SS = SolarSystemPMBH2(15 * AU, 0.5)   
    
    SS.addBody(Sun)
    SS.addBody(Earth)
    SS.addBody(Mercury)
    SS.addBody(Jupiter)
    SS.addBody(Mars)
    SS.addBody(Venus)
    
    
    
    for i in range(nbodies):
        belt = Planet(f'Belt {i}', 1., 1., 'str', 1., 1., 0.1, x_rand_belt[i], y_rand_belt[i], u_rand[i], v_rand[i])
        trojan = Planet(f'Trojan {i}', 1., 1., 'str', 1., 1., 0.1, x_rand_trojans[i], y_rand_trojans[i], u_rand[i], v_rand[i])
        SS.addBody(trojan)
        SS.addBody(belt)
        
        
    SS.Evolve(3600, method = 'Euler')
    pos = SS.getPositions()
    pos3.append(pos)
    

pos4 = []
    
for n in np.linspace(1, 3, 7, endpoint = True):
    
    rng1 = np.random.default_rng(seed = 82636473836)
    rng2 = np.random.default_rng(seed = 84634357238)
    
    Sun =       Planet('Sun',       sun['radius'],      1.3271244004193938E+11, 'red',      1.,          1.,        MSun,               0.,             0.,             0.,             0.)
    Earth =     Planet('Earth',     earth['radius'],    5.972e24 * G,           'blue',     1.,          1.,        earth['mass'],      earth['x'],     earth['y'],     earth['u'],     earth['v'])
    Jupiter =   Planet('Jupiter',   jupiter['radius'],  1.26686e8,              'red',      1.,          1.,        jupiter['mass'],    jupiter['x'],   jupiter['y'],   jupiter['u'],   jupiter['v'])
    Mercury =   Planet("Mercury",   2439.7,             0.3301e24 * G,          "gray",     57.909e6,    0.02056,   mercury['mass'],    mercury['x'],   mercury['y'],   mercury['u'],   mercury['v'])
    Mars =      Planet('Mars',      mars['radius'],     mars['mu'],             'green',    1.,          1.,        mars['mass'],       mars['x'],      mars['y'],      mars['u'],      mars['v'])
    Venus =     Planet('Venus',     venus['radius'],    venus['mu'],            'green',    1.,          1.,        venus['mass'],      venus['x'],     venus['y'],     venus['u'],     venus['v'])
    
    nbodies = int(10**n)
    rand_r_trojans = (rng1.normal(loc = 0., scale = 0.35, size = nbodies) + 5 ) * AU
    rand_theta_trojans = rng1.normal(loc = 0., scale = 0.5, size = nbodies) * np.pi / 6 + 0.63306 # math.atan2(Jupiter['y'], Jupiter['x'])     
    rand_theta_trojans[:nbodies // 2] += np.pi / 3
    rand_theta_trojans[nbodies // 2:9 * (nbodies // 10)] -= np.pi / 3
    rand_theta_trojans[9 * (nbodies // 10):] += np.pi
    
    rand_r_belt = (rng2.normal(loc = 0., scale = 0.5, size = nbodies ) + 2.66 ) * AU
    rand_theta_belt = rng2.random(nbodies) * 2 *  np.pi 
    
    x_rand_trojans = rand_r_trojans * np.cos(rand_theta_trojans)
    y_rand_trojans = rand_r_trojans * np.sin(rand_theta_trojans)
    
    x_rand_belt = rand_r_belt * np.cos(rand_theta_belt)
    y_rand_belt = rand_r_belt * np.sin(rand_theta_belt)
    
    speed_rand_trojans = np.sqrt(G * MSun / rand_r_trojans) 
    
    u_rand = -np.sin(rand_theta_trojans) * speed_rand_trojans
    v_rand = np.cos(rand_theta_trojans) * speed_rand_trojans
    
    speed_rand_belt = np.sqrt(G * MSun / rand_r_belt) 
    
    u_rand_belt = -np.sin(rand_theta_belt) * speed_rand_belt
    v_rand_belt = np.cos(rand_theta_belt) * speed_rand_belt
    
    mass_rand = 10**(np.random.rand(nbodies) * 18)
    
    SS = SolarSystemPMBH2(15 * AU, 1)   
    
    SS.addBody(Sun)
    SS.addBody(Earth)
    SS.addBody(Mercury)
    SS.addBody(Jupiter)
    SS.addBody(Mars)
    SS.addBody(Venus)
    
    
    
    for i in range(nbodies):
        belt = Planet(f'Belt {i}', 1., 1., 'str', 1., 1., 0.1, x_rand_belt[i], y_rand_belt[i], u_rand[i], v_rand[i])
        trojan = Planet(f'Trojan {i}', 1., 1., 'str', 1., 1., 0.1, x_rand_trojans[i], y_rand_trojans[i], u_rand[i], v_rand[i])
        SS.addBody(trojan)
        SS.addBody(belt)
        
        
    SS.Evolve(3600, method = 'Euler')
    pos = SS.getPositions()
    pos4.append(pos)
    
pos5 = []
    
for n in np.linspace(1, 3, 7, endpoint = True):
    
    rng1 = np.random.default_rng(seed = 82636473836)
    rng2 = np.random.default_rng(seed = 84634357238)
    
    Sun =       Planet('Sun',       sun['radius'],      1.3271244004193938E+11, 'red',      1.,          1.,        MSun,               0.,             0.,             0.,             0.)
    Earth =     Planet('Earth',     earth['radius'],    5.972e24 * G,           'blue',     1.,          1.,        earth['mass'],      earth['x'],     earth['y'],     earth['u'],     earth['v'])
    Jupiter =   Planet('Jupiter',   jupiter['radius'],  1.26686e8,              'red',      1.,          1.,        jupiter['mass'],    jupiter['x'],   jupiter['y'],   jupiter['u'],   jupiter['v'])
    Mercury =   Planet("Mercury",   2439.7,             0.3301e24 * G,          "gray",     57.909e6,    0.02056,   mercury['mass'],    mercury['x'],   mercury['y'],   mercury['u'],   mercury['v'])
    Mars =      Planet('Mars',      mars['radius'],     mars['mu'],             'green',    1.,          1.,        mars['mass'],       mars['x'],      mars['y'],      mars['u'],      mars['v'])
    Venus =     Planet('Venus',     venus['radius'],    venus['mu'],            'green',    1.,          1.,        venus['mass'],      venus['x'],     venus['y'],     venus['u'],     venus['v'])
    
    nbodies = int(10**n)
    rand_r_trojans = (rng1.normal(loc = 0., scale = 0.35, size = nbodies) + 5 ) * AU
    rand_theta_trojans = rng1.normal(loc = 0., scale = 0.5, size = nbodies) * np.pi / 6 + 0.63306 # math.atan2(Jupiter['y'], Jupiter['x'])     
    rand_theta_trojans[:nbodies // 2] += np.pi / 3
    rand_theta_trojans[nbodies // 2:9 * (nbodies // 10)] -= np.pi / 3
    rand_theta_trojans[9 * (nbodies // 10):] += np.pi
    
    rand_r_belt = (rng2.normal(loc = 0., scale = 0.5, size = nbodies ) + 2.66 ) * AU
    rand_theta_belt = rng2.random(nbodies) * 2 *  np.pi 
    
    x_rand_trojans = rand_r_trojans * np.cos(rand_theta_trojans)
    y_rand_trojans = rand_r_trojans * np.sin(rand_theta_trojans)
    
    x_rand_belt = rand_r_belt * np.cos(rand_theta_belt)
    y_rand_belt = rand_r_belt * np.sin(rand_theta_belt)
    
    speed_rand_trojans = np.sqrt(G * MSun / rand_r_trojans) 
    
    u_rand = -np.sin(rand_theta_trojans) * speed_rand_trojans
    v_rand = np.cos(rand_theta_trojans) * speed_rand_trojans
    
    speed_rand_belt = np.sqrt(G * MSun / rand_r_belt) 
    
    u_rand_belt = -np.sin(rand_theta_belt) * speed_rand_belt
    v_rand_belt = np.cos(rand_theta_belt) * speed_rand_belt
    
    mass_rand = 10**(np.random.rand(nbodies) * 18)
    
    SS = SolarSystemPMBH2(15 * AU, 3)   
    
    SS.addBody(Sun)
    SS.addBody(Earth)
    SS.addBody(Mercury)
    SS.addBody(Jupiter)
    SS.addBody(Mars)
    SS.addBody(Venus)
    
    
    
    for i in range(nbodies):
        belt = Planet(f'Belt {i}', 1., 1., 'str', 1., 1., 0.1, x_rand_belt[i], y_rand_belt[i], u_rand[i], v_rand[i])
        trojan = Planet(f'Trojan {i}', 1., 1., 'str', 1., 1., 0.1, x_rand_trojans[i], y_rand_trojans[i], u_rand[i], v_rand[i])
        SS.addBody(trojan)
        SS.addBody(belt)
        
        
    SS.Evolve(3600, method = 'Euler')
    pos = SS.getPositions()
    pos5.append(pos)
    
    
    
err2 = []
err3 = []
err4 = []
err5 = []

for i in range(7):
    err = np.sum(np.sqrt(np.sum((pos2[i] - pos1[i])**2, axis = 1))) / len(pos1[i])
    err2.append(err)
    err = np.sum(np.sqrt(np.sum((pos3[i] - pos1[i])**2, axis = 1))) / len(pos1[i])
    err3.append(err)
    err = np.sum(np.sqrt(np.sum((pos4[i] - pos1[i])**2, axis = 1))) / len(pos1[i])
    err4.append(err)
    err = np.sum(np.sqrt(np.sum((pos5[i] - pos1[i])**2, axis = 1))) / len(pos1[i])
    err5.append(err)
    


fig, ax = plt.subplots(figsize=(4.8, 4.8))
traj = np.loadtxt('run2.txt', delimiter=' ')
ax.scatter(traj[0, 0], traj[0, 1], color='yellow', marker='o', s = 50, label='Sun')
ax.scatter(traj[0, 4], traj[0, 5], color='brown', marker='o', s = 10, label='Mercury')
ax.scatter(traj[0,10], traj[0, 11], color='gray', marker='o', s = 20, label='Venus')
ax.scatter(traj[0, 2], traj[0, 3], color='green', marker='o', s = 20, label='Earth')
ax.scatter(traj[0, 8], traj[0, 9], color='red', marker='o', s = 15, label='Mars')
ax.scatter(traj[0, 6], traj[0, 7], color='purple', marker='o', s = 35, label='Jupiter')



ax.scatter(traj[0, 12], traj[0, 13], color='black', marker='.', s = 5, label='Asteroids')
for i in range(7,406):
    ax.scatter(traj[0, 2 * i], traj[0, 2*i + 1], color='black', marker='.', s = 5)

ax.legend(loc = 'upper left')
ax.set_xlim(-7.5*AU, 7.5*AU)
ax.set_ylim(-7.5*AU, 7.5*AU)
plt.show()

        

   ''' 

'''
    SS = SolarSystemPMBH2(15 * AU, theta)
    
    SS.addBody(Sun)
    SS.addBody(Earth)
    SS.addBody(Mercury)
    SS.addBody(Jupiter)
    SS.addBody(Mars)
    SS.addBody(Venus)
    
    # Testing Quadtree implementation
    #SS.buildTree()
    
    #displaySystem(SS, 1)
    
    
    # Testing force computation
    
    #SS.Orbit(12 * 365, 12 * 365*24*3600, method ='Euler')
    rng1 = np.random.default_rng(seed = 82636473836)
    rng2 = np.random.default_rng(seed = 84634357238)
    
    nbodies = 400
    rand_r_trojans = (rng1.normal(loc = 0., scale = 0.35, size = nbodies) + 5 ) * AU
    rand_theta_trojans = rng1.normal(loc = 0., scale = 0.5, size = nbodies) * np.pi / 6 + 0.63306 # math.atan2(Jupiter['y'], Jupiter['x'])     
    rand_theta_trojans[:nbodies // 2] += np.pi / 3
    rand_theta_trojans[nbodies // 2:9 * (nbodies // 10)] -= np.pi / 3
    rand_theta_trojans[9 * (nbodies // 10):] += np.pi
    
    rand_r_belt = (rng2.normal(loc = 0., scale = 0.4, size = nbodies ) + 2.66 ) * AU
    rand_theta_belt = rng2.random(nbodies) * 2 *  np.pi 
    
    x_rand_trojans = rand_r_trojans * np.cos(rand_theta_trojans)
    y_rand_trojans = rand_r_trojans * np.sin(rand_theta_trojans)
    
    x_rand_belt = rand_r_belt * np.cos(rand_theta_belt)
    y_rand_belt = rand_r_belt * np.sin(rand_theta_belt)
    
    speed_rand_trojans = np.sqrt(G * MSun / rand_r_trojans) 
    
    u_rand = -np.sin(rand_theta_trojans) * speed_rand_trojans
    v_rand = np.cos(rand_theta_trojans) * speed_rand_trojans
    
    speed_rand_belt = np.sqrt(G * MSun / rand_r_belt) 
    
    u_rand_belt = -np.sin(rand_theta_belt) * speed_rand_belt
    v_rand_belt = np.cos(rand_theta_belt) * speed_rand_belt
    
    mass_rand = 10**(np.random.rand(nbodies) * 18)
    
    for i in range(nbodies):
        #belt = Body(1e18, x_rand_belt[i], y_rand_belt[i], u_rand[i], v_rand[i])
        trojan = Body(1e18, x_rand_trojans[i], y_rand_trojans[i], u_rand[i], v_rand[i])
        SS.addBody(trojan)
        #SS.addBody(belt)
        
    #SS.buildTree()
    #displaySystem(SS)
    #SS.Orbit(400 * 365 // 20,  400 * 365 * 24 * 3600, method = 'Leapfrog')
    
    # theta vs runtime
    
    t1 = time.time()
    SS.Evolve(3600)
    t2 = time.time()
    t_theta1.append(t2 - t1)
'''
'''
run3 = np.zeros((400*365 // 20, 406, 2))
for i, body in enumerate(SS.getBodies()):
    run3[:, i] = body.trajectory
run3 = run3.reshape(run3.shape[0], -1)
np.savetxt('run3.txt', run3)
'''
'''
t_PP = []
t_BH = []
# Testing Large Scale Particles
for n in np.linspace(1, 4, 10, endpoint = True):
    nbodies = int(10**n)
    rand_r_trojans = (np.random.normal(loc = 0., scale = 0.35, size = nbodies) + 5 ) * AU
    rand_theta_trojans = np.random.normal(loc = 0., scale = 0.5, size = nbodies) * np.pi / 6 + 0.63306 # math.atan2(Jupiter['y'], Jupiter['x'])     
    rand_theta_trojans[:nbodies // 2] += np.pi / 3
    rand_theta_trojans[nbodies // 2:9 * (nbodies // 10)] -= np.pi / 3
    rand_theta_trojans[9 * (nbodies // 10):] += np.pi
    
    rand_r_belt = (np.random.normal(loc = 0., scale = 0.5, size = nbodies ) + 2.66 ) * AU
    rand_theta_belt = np.random.rand(nbodies) * 2 *  np.pi 
    
    x_rand_trojans = rand_r_trojans * np.cos(rand_theta_trojans)
    y_rand_trojans = rand_r_trojans * np.sin(rand_theta_trojans)
    
    x_rand_belt = rand_r_belt * np.cos(rand_theta_belt)
    y_rand_belt = rand_r_belt * np.sin(rand_theta_belt)
    
    speed_rand_trojans = np.sqrt(G * MSun / rand_r_trojans) 
    
    u_rand = -np.sin(rand_theta_trojans) * speed_rand_trojans
    v_rand = np.cos(rand_theta_trojans) * speed_rand_trojans
    
    speed_rand_belt = np.sqrt(G * MSun / rand_r_belt) 
    
    u_rand_belt = -np.sin(rand_theta_belt) * speed_rand_belt
    v_rand_belt = np.cos(rand_theta_belt) * speed_rand_belt
    
    mass_rand = 10**(np.random.rand(nbodies) * 18)
    
    SS1 = SolarSystemBH(15 * AU, 0.5)
    SS2 = SolarSystemPP(15 * AU)

    SS1.addBody(Sun)
    SS1.addBody(Earth)
    SS1.addBody(Mercury)
    SS1.addBody(Jupiter)
    SS1.addBody(Mars)
    SS1.addBody(Venus)
    
    SS2.addBody(Sun)
    SS2.addBody(Earth)
    SS2.addBody(Mercury)
    SS2.addBody(Jupiter)
    SS2.addBody(Mars)
    SS2.addBody(Venus)
    
    for i in range(nbodies):
        belt = Planet(f'Belt {i}', 1., 1., 'str', 1., 1., 0.1, x_rand_belt[i], y_rand_belt[i], u_rand[i], v_rand[i])
        trojan = Planet(f'Trojan {i}', 1., 1., 'str', 1., 1., 0.1, x_rand_trojans[i], y_rand_trojans[i], u_rand[i], v_rand[i])
        SS1.addBody(trojan)
        SS1.addBody(belt)
        SS2.addBody(trojan)
        SS2.addBody(belt)
        
        
    t1 = time.time()
    SS1.Evolve(3600, 0.5)
    t2 = time.time()
    SS2.Evolve(3600)
    t3 = time.time()
    
    t_BH.append(t3 - t2)
    t_PP.append(t2 - t1)
    print(t3 - t2, t2 - t1)
    
    

#SS.buildTree()
#displaySystem(SS)

# Test times

    '''
'''
dom = np.linspace(10, 10 ** 3.5, 100)

plt.figure()
plt.loglog(2 * 10**(np.linspace(1, 3.33, 8, endpoint = True)), t_BH, 'r.', label = 'Barnes-Hut')
plt.loglog(2 * 10**(np.linspace(1, 3.33, 8, endpoint = True)), t_PP, 'b.', label = 'Particle-Particle')
plt.loglog(dom, dom**2 / 10000, 'r--', label = '$\sim N^2$')
plt.loglog(dom, dom * np.log10(dom) / 3000, 'b--', label = '$\sim NlogN$')
plt.xlabel('$N$')
plt.ylabel('$t_{computation}$')
#plt.title('Asymptotic Runtime: Barnes-Hut vs. Particle-Particle Method')
plt.legend()
plt.show()
'''
'''
dom = np.linspace(10, 10 ** 3.5, 100)

plt.figure()
plt.loglog(2 * 10**(np.linspace(1, 3, 7, endpoint = True)), t_theta1, '.', c = 'C1', label = r'$\theta = 0$')
plt.loglog(2 * 10**(np.linspace(1, 3, 7, endpoint = True)), t_theta2, 'b.', label = r'$\theta = 0.5$')
plt.loglog(2 * 10**(np.linspace(1, 3, 7, endpoint = True)), t_theta3, 'g.', label = r'$\theta = 1$')
plt.loglog(2 * 10**(np.linspace(1, 3, 7, endpoint = True)), t_theta4, 'm.', label = r'$\theta = 3$')
plt.loglog(dom, dom**2 / 5000, '--', c='C1', label = '$\propto N^2$')
plt.loglog(dom, dom / 3000, 'm--', label = '$\propto N$')
plt.loglog(dom, dom * np.log10(dom) / 1500, 'b--', label = '$\propto N\log N$')
plt.xlabel('$N$')
plt.ylabel('$t_{computation}$')
#plt.title('Asymptotic Runtime: Barnes-Hut vs. Particle-Particle Method')
plt.legend()
plt.show()
'''
'''
dom = np.linspace(0.2, 3, 100)

plt.figure()
plt.plot(10 **np.linspace(1, 3, 7, endpoint = True), err2, 'c.', label = r'$\theta = 0$')
plt.plot(10 **np.linspace(1, 3, 7, endpoint = True), err3, 'm.', label = r'$\theta = 0.5$')
plt.plot(10 **np.linspace(1, 3, 7, endpoint = True), err4, 'g.', label = r'$\theta = 1$')
plt.plot(10 **np.linspace(1, 3, 7, endpoint = True), err5, 'b.', label = r'$\theta = 3$')
#plt.plot(dom, np.exp(-dom)
plt.xlabel('Number of bodies $N$ ')
plt.ylabel('Mean Error relative to PP Method')
#plt.title('Asymptotic Runtime: Barnes-Hut vs. Particle-Particle Method')
plt.legend()
plt.show()


#SS.Orbit(12 * 365, 12 * 365.25 * 24 * 3600, method = 'Leapfrog')
'''
'''
plt.figure()
for body in SS.getBodies():
    if np.sqrt(body.trajectory[-1, 0]**2 +  body.trajectory[-1, 1]**2) < 20 * AU:
        plt.plot(body.trajectory[:, 0], body.trajectory[:, 1])
  


#plt.axes('equal')    
plt.show()
'''

#SS.Evolve(24*3600, method = 'Euler')

#SS.Orbit(365, 365*24*3600, method ='Euler')


#SS.Orbit(730, 2 * 365 * 24 * 3600)
'''

nbodies = 100
pos_rand = np.random.rand(2, nbodies) - 0.5
pos_rand *=  AU 
for i in range(nbodies):
    body = Planet(f'{i}', 1., 1., 'str', 1., 1., 1e15, pos_rand[0, i], pos_rand[1, i], 0., 0.)
    SS.addBody(body)
    
SS.buildTree()
displaySystem(SS)
SS.Orbit(365, 365 * 24 * 3600)
'''
#for body in SS.getBodies():
   # SS.addNode(body)
#SS.Orbit(100000, 5000000000, 'Euler', GR = True)

# E_x = Earth.trajectory[:,0]
# E_y = Earth.trajectory[:,1]

#M_x = Mercury.trajectory[:,0]
#M_y = Mercury.trajectory[:,1]



# SS.initGrids()
# SS.gridParticles()




# # 3D end state behaviours
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.plot(E_x, E_y, np.zeros(len(Earth.trajectory[:,1])))
# ax.plot(J_x, J_y,  np.zeros(len(Earth.trajectory[:,1])))
# ax.scatter(0.,0.,0., s= 300, color = 'red')

# plt.show()


# Testing times




# 2D animation of evolution
import matplotlib.animation
# The animation requires FFmpeg 

#fig, ax = plt.subplots()

def animate(t):
    plt.cla()
    plt.xlim(-0.5*AU, 0.5*AU)
    plt.ylim(-0.5*AU, 0.5*AU)

    # plt.plot(E_x[0: t*500], E_y[0: t*500])
    #plt.plot(M_x[0: t*500], M_y[0: t*500])

    plt.scatter(0.,0., s= 200, color = 'red')



#plt.rcParams["animation.html"] = "html5"
 
#anim = matplotlib.animation.FuncAnimation(fig, animate, frames=199, interval = 80)


#anim

#anim.save("test.mp4")
