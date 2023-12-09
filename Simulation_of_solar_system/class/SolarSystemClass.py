# -*- coding: utf-8 -*-
from PlanetClass import Body, Planet
import numpy as np

class SolarSystem:
    
    def __init__(self, size):
        self.size = size   # cartesian dimesions of 2D plane
        self.bodies = []   # stores all N-bodies
        self.G = 6.67e-20  # gravitational constant [km^3 kg^-1 s^-2]
        self.ep = 1e-4     # epsilon (softening parameter)
        
    def addBody(self, body: Body):
        self.bodies.append(body)
        
    def getBodies(self):
        return self.bodies
    
    def getMasses(self):
        return np.array(list(map(lambda b: b.getMass(), self.bodies)))
        
    def getVelocities(self):
        return np.array(list(map(lambda b: b.getVel(), self.bodies)))
    
    def getPositions(self):
        return np.array(list(map(lambda b: b.getPos(), self.bodies)))
    
    def updatePositions(self, pos):
        for b, xy in zip(self.bodies, pos):
            b.setPos(xy)
        
    def updateVelocities(self, vels):
        for b, uv in zip(self.bodies, vels):
            b.setVel(uv)
            
    def initTrajectories(self, nsteps, dim = 2):
        for b in self.bodies:
            b.initTrajectory(nsteps, dim)
            
    def storeTrajectories(self, step):
        for b in self.bodies:
            b.storeTrajectory(step)

    # this part is for RK4 only
    def init_tras(self):
        for b in self.bodies:
            b.init_tra_vary()
    def append_tras(self):
        for b in self.bodies:
            b.append_tra()
    def reverse_tras(self):
        for b in self.bodies:
            b.reverse_tra()

        

def getDis(x1, x2):
    r = x2 - x1
    rr = np.linalg.norm(r)
    r_hat = r / rr    
    return rr, r_hat


def get_accs(self, dt, velo):
    bodies = self.bodies

    pos_list = self.getPositions() + dt * velo
    
    accs = np.zeros((len(bodies), 2))
    
    for i, curBody in enumerate(bodies):

        a = np.array([0., 0.])

        for k, body in enumerate(bodies):
            if body != curBody:
                rr, r_hat = getDis(pos_list[i], pos_list[k])
                a += body.getMu() * r_hat / (rr**2 + self.ep**2)
                # softened force ep = 1e-4

        accs[i] = a      

    return accs     

def get_accs_GR (self, dt, velo):
    bodies = self.bodies

    pos_list = self.getPositions() + dt * velo
    
    accs = np.zeros((len(bodies), 2))
    
    for i, curBody in enumerate(bodies):

        a = np.array([0., 0.])

        for k, body in enumerate(bodies):
            if body != curBody:
                rr, r_hat = getDis(pos_list[i], pos_list[k])
                a += body.getMu() * r_hat / (rr**2 + self.ep**2) * (1 - 2/(9e10 * rr))
                # softened force ep = 1e-4

        accs[i] = a      

    return accs   



class SolarSystemPP(SolarSystem):
    """ Particle-Particle (naive) simulation"""
        
    def __init__(self, size):
            super().__init__(size)
    
    # should make evolve a method of SolarSystem accessible to all subclasses, maybe put accs as arg

    def Evolve(self, dt, method = 'Euler'):
        """ Perform an integration step (particle-particle) on the N-body system according to specified integrator method
        """
        
        v_0 = self.getVelocities()

        accs_0 = get_accs(self, dt = 0, velo = v_0)
            
        if method == 'Euler':
            # Perform Euler integration step
            newVels = v_0 + accs_0 * dt
            self.updateVelocities(newVels)
            
            newPos = self.getPositions() + newVels * dt
            self.updatePositions(newPos)

        if method == 'RK4':
            # Perform 4th order Runge-Kutta integration step
            v_1 = v_0 + accs_0 * dt/2
            a_1 = get_accs(self, dt = dt/2, velo = v_1)

            v_2 = v_0 + a_1 * dt/2
            a_2 = get_accs(self, dt = dt/2, velo = v_2)

            v_3 = v_0 + a_2 * dt
            a_3 = get_accs(self, dt = dt, velo = v_3)

            v_4 = v_0 + dt/6 * (accs_0 + 2*a_1 + 2*a_2 + a_3)

            self.updateVelocities(v_4)
            
            newPos = self.getPositions() + v_4 * dt
            self.updatePositions(newPos)

    def Evolve_GR(self, dt, method = 'Euler'):
        """ Perform an integration step (particle-particle) on the N-body system according to specified integrator method
            while accounting for general relativity correction factor
        """

        v_0 = self.getVelocities()

        accs_0 = get_accs_GR(self, dt = 0, velo = v_0)
            
        if method == 'Euler':
            # Perform Euler integration step
            newVels = v_0 + accs_0 * dt
            self.updateVelocities(newVels)
            
            newPos = self.getPositions() + newVels * dt
            self.updatePositions(newPos)

        if method == 'RK4':
            # Perform 4th order Runge-Kutta integration step
            v_1 = v_0 + accs_0 * dt/2
            a_1 = get_accs_GR(self, dt = dt/2, velo = v_1)

            v_2 = v_0 + a_1 * dt/2
            a_2 = get_accs_GR(self, dt = dt/2, velo = v_2)

            v_3 = v_0 + a_2 * dt
            a_3 = get_accs_GR(self, dt = dt, velo = v_3)

            v_4 = v_0 + dt/6 * (accs_0 + 2*a_1 + 2*a_2 + a_3)

            self.updateVelocities(v_4)
            
            newPos = self.getPositions() + v_4 * dt
            self.updatePositions(newPos)
            
            
    def Orbit(self, nsteps, t_end, method, t0 = 0, GR = False):
        """ Perform integration on the N-body system for the entire period according to specified integrator method
        """
        dt = (t_end - t0) / nsteps
        self.initTrajectories(nsteps)
        
        self.storeTrajectories(0)
        
        if GR == False:
            for i in range(1, nsteps):
                self.Evolve(dt, method)
                self.storeTrajectories(i)
        else:
            for i in range(1, nsteps):
                self.Evolve_GR(dt, method)
                self.storeTrajectories(i)




    """
    The following is a specical case that the function 
    are tailored to perform adaptive RK4 methods with GR
    """

    def Evolve_RK_GR(self, dt):
        v_0 = self.getVelocities()
        accs_0 = get_accs_GR(self, dt = 0, velo = v_0)
        v_1 = v_0 + accs_0 * dt/2
        a_1 = get_accs_GR(self, dt = dt/2, velo = v_1)
        v_2 = v_0 + a_1 * dt/2
        a_2 = get_accs_GR(self, dt = dt/2, velo = v_2)
        v_3 = v_0 + a_2 * dt
        a_3 = get_accs_GR(self, dt = dt, velo = v_3)
        v_4 = v_0 + dt/6 * (accs_0 + 2*a_1 + 2*a_2 + a_3)
        newPos = self.getPositions() + v_4 * dt
        self.updateVelocities(v_4)
        self.updatePositions(newPos)
        return newPos



    def Orbit_RK_GR(self, t_end, tol = 10):
        self.init_tras()

        t_end = t_end 
        tol = tol
        t = 0.0
        dt = 1000
        
        while t < t_end:

            while True:
                x00 = self.getPositions()
                v00 = self.getVelocities()

                x1 = self.Evolve_RK_GR(dt)
                self.updatePositions(x00)
                self.updateVelocities(v00)

                x2 = self.Evolve_RK_GR(dt/2)
                x3 = self.Evolve_RK_GR(dt/2)
                err = np.max(abs((x3-x1)))

                if err < tol:
                    t = t + dt
                    dt = dt * 2
                    # check to see whether this step size would
                    # take us beyond t=t_end and if so adjust accordingly
                    if t + dt > t_end:
                        dt = t_end - t
                    # store 2-step evolution
                    self.append_tras()
                    break
                else:
                    dt = dt / 2
                    self.updatePositions(x00)
                    self.updateVelocities(v00)
        
        self.reverse_tras()



        


