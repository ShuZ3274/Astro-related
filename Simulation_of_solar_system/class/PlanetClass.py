# -*- coding: utf-8 -*-
import numpy as np


class Body:
    
    def __init__(self, mass: float, x: float, y: float, u: float, v: float):
        self.m = mass # [kg]
        self.x = x    # [km]
        self.y = y    # [km]
        self.u = u    # [km/s]
        self.v = v    # [km/s]
        self.trajectory = None # stores orbit trajectory (used in SolarSystem classes)
        
    def __str__(self):
        return f'Position: ({self.x}, {self.y})' + '\n' + f'Velocity: ({self.u}, {self.v})' + '\n' + f'Mass: {self.m}'
        
    def getMass(self):
        return self.m
         
    def getPos(self):
        return np.array([self.x, self.y])
    
    def setPos(self, pos):
        x, y = pos[0],  pos[1]
        self.setX(x)
        self.setY(y)
        
    def setX(self, x):
        self.x = x
        
    def setY(self, y):
        self.y = y
        
    def getVel(self):
        return np.array([self.u, self.v])
    
    def setVel(self, vels):
        u, v = vels[0],  vels[1]
        self.setU(u)
        self.setV(v)
        
    def setU(self, u):
        self.u = u
        
    def setV(self, v):
        self.v = v
        
    def getDistance(self, body):
        x1 = self.getPos()
        x2 = body.getPos()
        
        r = x2 - x1
        rr = np.linalg.norm(r) # euclidian distance
        r_hat = r / rr # unit vector
        return rr, r_hat
    
    def initTrajectory(self, nsteps, dim = 2):
        self.trajectory = np.zeros((nsteps, dim))
        
    def storeTrajectory(self, step):
        self.trajectory[step] = self.getPos()
    

    # the following functions are the tailored for RK4
    def init_tra_vary(self, dim = 2):
        self.trajectory = np.array([self.getPos()])
    def append_tra(self):
        self.trajectory = np.insert(self.trajectory, 0, self.getPos(), axis=0)
    def reverse_tra(self):
        self.trajectory = self.trajectory[::-1]



class Planet(Body):
        
    def __init__(self, name: str, radius: float, mu: float, traj_color: str, 
                 semi_maj_axis: float, eccentricity: float, mass, x, y, u, v):
        super().__init__(mass, x, y, u, v)


        self.name = name
        self.radius = radius
        self.traj_c = traj_color
        self.mu = mu # GM
        self.maj_ax = semi_maj_axis # known orbit
        self.e = eccentricity       # known orbit
        
    def __str__(self):
        return f'Name: {self.name}' + '\n' +f'Radius: {self.radius}' + '\n'+ super().__str__() 
    
    def getMu(self):
        return self.mu
    
    def getRadius(self):
        return self.radius

        


