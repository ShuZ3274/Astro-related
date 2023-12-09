# -*- coding: utf-8 -*-
from PlanetClass import Body, Planet
from SolarSystemClass import SolarSystem, getDis
import numpy as np


class SolarSystemBH(SolarSystem):
    """ Barnes-Hut simulation"""
    
    def __init__(self, size, theta, G = 6.67e-20, EP = 1e-4):
        super().__init__(size)
        self.root = Node(bbox = [-size / 2, size / 2, -size / 2, size / 2]) # dimensions of 2D plane [xmin, xmax, ymin, ymax]
        self.theta = theta           # ratio used for CM approximation
        self.G = G                   # gravitational constant [km^3 kg^-1 s^-2]
        self.EP = EP                 # epsilon (softening parameter)
        
    def initTree(self):
        size = self.size
        self.root = Node(bbox = [-size / 2, size / 2, -size / 2, size / 2])
        
        
    def buildTree(self):
        for body in self.getBodies():
            x, y = body.getPos()
            mass = body.getMass()
            xmin, xmax, ymin, ymax = self.root.bbox
            if x < xmin or x > xmax or y < ymin or y > ymax: 
                return
            self.root.insert(x, y, mass)
        
        
        
    def getAcc(self, node, x, y, mass):
        # shouldnt ever happen but just in case
        if node.tot_mass is None:
            return np.array([0., 0.])
        # to remove self interactions
        if node.CM[0] == x and node.CM[1] == y and node.tot_mass == mass:
            return np.array([0., 0.])
        gridsize = node.bbox[1] - node.bbox[0]
        r, r_hat = getDis(np.array([x,y]), node.CM)
        # if ratio of gridsize of current node to distance from the body to the CM of the node
        # is smaller than theta (i.e. the group of particles in the node that we want to approximate as 
        # one has a relatively small diameter compared to the distance between the body and the 
        # group of particles) we go ahead with the approximation, alternatively, if the node were in has
        # only a single body in it (no children) we proceed with the same method but this calculation 
        # will now be equivalent to a particle-particle calculation.
        if gridsize / r < self.theta or node.children is None:
            return self.G * node.tot_mass / (r**2 + self.EP) * r_hat
        # if not go into each individual quadrant of that node (for better accuracy) and the  
        # total contribution will be the sum of the contributions from each quadrant, call getAcc() 
        # recursively on each quadrant (child), corresponds to a depth first search
        else:
            a_tot = np.array([0., 0.])
            for child in node.children:
                if child is not None:
                    a_tot += self.getAcc(child, x, y, mass)
        return a_tot
    
    def Evolve(self, dt, method = 'Leapfrog', step = 1):
        """ Perform an integration step (Barnes-Hut) on the N-body system according to specified integrator method
        """
        self.initTree()
        self.buildTree()
        
        # Calculate the accelerations on each body
        accs = np.zeros_like(self.getPositions())
        for i, body in enumerate(self.getBodies()):
            x, y = body.getPos()
            mass = body.getMass()
            accs[i] += self.getAcc(self.root, x, y, mass)
        
        pos = self.getPositions()
        vels = self.getVelocities()
        
        
        if method == 'Euler':
            # Perform Euler integration step
            newVels = vels + accs * dt
            self.updateVelocities(newVels)
            
            newPos = pos + newVels * dt
            self.updatePositions(newPos)

        '''if method == 'RK4':
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
            '''
            
        if method == 'Leapfrog':
            # Perform Leapfrog integration in two steps using Evolve recursively for convenience
           if step == 1:
                v_half = vels + accs * dt / 2
                newPos = pos + v_half * dt
                self.updatePositions(newPos)
                self.updateVelocities(v_half)
                self.Evolve(dt, method = 'Leapfrog', step = 2)
           elif step == 2:
               newVels = vels + accs * dt / 2
               self.updateVelocities(newVels)
               
    def Orbit(self, nsteps, t_end, method = 'Leapfrog', t0 = 0):
        """ Perform integration on the N-body system for the entire period according to specified integrator method
        """
        dt = (t_end - t0) / nsteps # timestep
        
        self.initTrajectories(nsteps)
        self.storeTrajectories(0) # initial position
        
        for i in range(1, nsteps):
            self.Evolve(dt, method)
            self.storeTrajectories(i)
            print(f'{i * 20 /365} years elapsed')
    
        
    def displayTree(self, node):
        if node is None:
            return
        elif node.children is None:
            print(node)
            return
        else:
            for child in node.children:
                self.displayTree(child)
        
    

class Node:
    
    def __init__(self, tot_mass = None, CM = None, bbox = None):
        self.tot_mass = tot_mass # cumulative mass of subtree branching from self
        self.CM = CM             # center of mass 
        self.bbox = bbox         # dimensions of quadrant (cartesian)
        self.children = None     # None if external/leaf node
        
    def __str__(self):
        return f'Total Mass: {self.tot_mass} \n Center of Mass: {self.CM} \n Box: {self.bbox}'
        
    def initChildren(self):
        """ Creates quad children """
        self.children = [None, None, None, None]
    
    
    def insert(self, x, y, mass):
        """ Recursively inserts new body into quad tree branching from self 
            x,y,m : parameters (2d cartesian coordinates and mass) of new body
        """
        
        # if the node is empty, give it the attributes of the new body
        if self.tot_mass is None:
            self.tot_mass = mass
            self.CM = np.array([x, y])
            return # terminate recursion, we have reached a leaf
        
        # if the node is external/leaf (and already contains a body), 
        # subdivide region and recursively insert new particle
        elif self.children is None:

            self.initChildren() 
            
            # move the body originally occupying the node to the appropriate quadrant/child
            cur_quad = self.getQuad(*self.CM)
            cur_quad_bbox = self.constructBBOX(cur_quad)
            self.children[cur_quad] = Node(bbox = cur_quad_bbox)
            self.children[cur_quad].insert(*self.CM, self.tot_mass)
            
            # add the new body to the tree by adding it to the appropriate quadrant/child
            new_quad = self.getQuad(x, y)
            if self.children[new_quad] is None:
                new_quad_box = self.constructBBOX(new_quad)
                self.children[new_quad] = Node(bbox = new_quad_box)
            self.children[new_quad].insert(x, y, mass)
            
        # if node is internal/branch, recursively insert new body
        else:
            new_quad = self.getQuad(x, y)
            # check if there are other bodies in the same quadrant
            if self.children[new_quad] is None:
                # if not, create quadrant
                new_quad_bbox = self.constructBBOX(new_quad)
                self.children[new_quad] = Node(bbox = new_quad_bbox)
            # add body to quadrant
            self.children[new_quad].insert(x, y, mass)
            
        # update info stored in node (center of mass and total mass) to account for the new body
        xcm = self.CM[0] * self.tot_mass + x * mass
        ycm = self.CM[1] * self.tot_mass + y * mass
        self.tot_mass += mass
        xcm /= self.tot_mass
        ycm /= self.tot_mass 
        self.CM[0] = xcm
        self.CM[1] = ycm
            
           
    def getQuad(self, x, y):
        """ Get node's quadrant/child to which (x,y) coordinates belongs 
            (index {0,1,2,3} uses cartesian plane convention)
        """
        xmin, xmax, ymin, ymax = self.bbox
        
        if x < xmin or x > xmax or y < ymin or y > ymax: # if outside box
            raise IndexError(f'Position exceeds bbox of current node, \n Position = ({x},{y}), BBOX = ({self.bbox})')
            
        if y > (ymin + ymax) / 2: 
            if x > (xmin + xmax) / 2: 
                return 0 # upper right
            else: 
                return 1 # upper left
        elif x > (xmin + xmax) /2: 
            return 3     # lower right
        else: 
            return 2     # lower left

        
    def constructBBOX(self, quad):
        """ Calculate (cartesian) bounds of new subquadrant 
        """
        xmin, xmax, ymin, ymax = self.bbox
        if quad == 0:
            return [(xmin + xmax) /2, xmax, (ymin + ymax) /2, ymax]
        elif quad == 1:
            return [xmin, (xmin + xmax) /2, (ymin + ymax) /2, ymax]
        elif quad == 2:
            return [xmin, (xmin + xmax) /2, ymin, (ymin + ymax) /2]
        elif quad == 3:
            return [(xmin + xmax) /2, xmax, ymin, (ymin + ymax) /2]
        else: 
            raise IndexError('Quadrant indices must be int between 0-3')
            
     
            
# Visualization
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def plotQuadtree(ax, node):
    if node is None:
        return
    elif node.children is None:
        # leaf -> plot bbox
        bbox = node.bbox
        rect = Rectangle((bbox[0], bbox[2]), bbox[1] - bbox[0], bbox[3] - bbox[2], fill=False, edgecolor='blue')
        ax.add_patch(rect)
        return
    else:
        # internal node -> recursively plot children
        for child in node.children:
            plotQuadtree(ax, child)

def plotParticles(ax, system):
    positions = system.getPositions()
    ax.scatter(positions[6:, 0], positions[6:, 1], color='black', marker='.', s = 5, label='Asteroids')

def plotCelestials(ax, system):
    positions = system.getPositions()
    ax.scatter(positions[0, 0], positions[0, 1], color='yellow', marker='o', s = 50, label='Sun')
    ax.scatter(positions[2, 0], positions[2, 1], color='brown', marker='o', s = 10, label='Mercury')
    ax.scatter(positions[5, 0], positions[5, 1], color='gray', marker='o', s = 20, label='Venus')
    ax.scatter(positions[1, 0], positions[1, 1], color='green', marker='o', s = 20, label='Earth')
    ax.scatter(positions[4, 0], positions[4, 1], color='red', marker='o', s = 15, label='Mars')
    ax.scatter(positions[3, 0], positions[3, 1], color='purple', marker='o', s = 35, label='Jupiter')
  
    
def displaySystem(system):
    fig, ax = plt.subplots(figsize=(6.4, 6.4))
    ax.set_xlim(system.root.bbox[0], system.root.bbox[1])
    ax.set_ylim(system.root.bbox[2], system.root.bbox[3])
    ax.set_xlabel('x[km]')
    ax.set_ylabel('y[km]')
    
    # plot particles
    plotParticles(ax, system)
    # plot segmentation
    plotQuadtree(ax, system.root)
    
    plotCelestials(ax, system)

    ax.legend()
    plt.show()



            
            
        
        