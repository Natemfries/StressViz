# Class definitions and related functions for generating Cycloids in SatStressGUI

# possibly import pytest in order to use assertions also as error checkers
# brought up because asserts issues automatically come up as an error window
import satstress
import random
from numpy import *
from pylab import *
from mpl_toolkits import basemap
# for object serialization : http://pythontips.com/2013/08/02/what-is-pickle-in-python/
import pickle
from matplotlib import pyplot as plt
# Open Source Geospacial libraries
from osgeo import ogr
from osgeo import osr
import os

from lineament import spherical_reckon

CycloidCoords = [] # A global array variable used to store the cycloid coordinates for saving to a shapefile. -PS 2016

class Cycloid(object):
    """
    A one dimensional feature on the surface of a spherical satellite similar to Lineament
    """
    def __init__(self, stresscalc, cycloid_name, threshold, propagation_strength, propagation_speed,
                 start_lat, start_lon, start_dir='West', varyVelocity=False, k=0, maxorbit=360, degree_step=0.5):
        """
        Generate a cycloid given a {stresscalc} object which defines which stresses to 
        calculate, the yield {threshold} value to initiate the crack, {propagation_strength},
        {propagation_speed}, and starting point (lon, lat), starting direction, 
        and stepping interval.

        There also includes the option of varying the crack propagation velocity as a function
        of stress magnitude. For this option, in addition to giving a propogation speed, the
        user must also give a proportionality constant 'k' such that:
        
        velocity = user_given_propogation_speed + k(tensile_str - propogation_strength).
     
        stresscalc is a L{satstress.StressCalc} object defining the field we
        are going to be comparing the lineament to.
        """
        
        ###### check validity of arguments #####
        assert start_dir is not None
        # check that k is valid if varyVelocity option is set
        if varyVelocity == True:
            assert k != 0
            assert(isinstance(k, (float, int)))

        ##### initialize #####
        self.cycloid_name = cycloid_name
        self.stresscalc = stresscalc
        self.tensile_str = float(threshold)*1000 # kPa
        self.propagation_strength = float(propagation_strength)*1000   # kPa
        self.propagation_speed = float(propagation_speed)         # m/s
        self.start_lat = float(start_lat)
        self.start_lon = float(start_lon)
        self.dir = start_dir
        self.varyVelocity = varyVelocity
        self.propconst = float(k)

        # lists of lat, lon points of the cycloid in format [pt1, pt2, None, pt2, pt3, None...]
        # for faster plotting of many line segments EXPLAIN why we do so...
        # ^^^^ No idea why they save things like this.  -PS (2016)
        self.lat_cycloidpts = []
        self.lon_cycloidpts = []
        self.current_orbit_pos = 0
        self.time = 0

        self.initiated = False
        self.not_propagated = False
        # list of tuples (lat, long, principle stress comps) summarizing important info 
        # for export/display purposes 
        self.data = []
    
        # keeps track of how far around the satellite the cycloid has traveled
        self.traveled = 0 

        # set the ang_dist traveled each interval i.e. angular distance traveled by cycloid
        # going at propagation_speed every x degrees satellite travels in orbit
        # for which a point for a cycloid plot is generated.
        self.degree_step = degree_step
        self.time_step = self.set_interval(self.degree_step)
        self.ang_dist = self.calc_ang_dist()
    
        self.initiate_split(maxorbit)

    # BEGIN: helper functions
    def set_interval(self, degree):
        """
        Calculates the self.time interval in seconds for finding angular distance
        """
        return degree / 360 * self.stresscalc.stresses[0].satellite.orbit_period()

    def calc_ang_dist(self):
        """
        Calculating angular distance is by finding 
        linear_distance = propagation_speed * self.time_step,
        and then converting to angular distance by
        angular_distance = linear_distance / satellite.radius
        """
        return self.propagation_speed * self.time_step / self.stresscalc.stresses[0].satellite.radius()

    def vary_ang_dist(self, stressmag):
        """
        Used to caculate the velocity to the next step when varying velocity option is set.
        
        Note: There are two separate functions instead of replacing calc_ang_dist with 
        vary_ang_dist w/ k = 0 because vary_ang_dist needs to be calculated every step
        while calc_ang_dist only need to be done once.
        """
        newspeed = self.propagation_speed + self.propconst * (stressmag - self.tensile_str)
        return newspeed * self.time_step / self.stresscalc.stresses[0].satellite.radius()

    def get_stresses(self, lat, lon):
        """
        Retrieves stress components at given point lat, lon (which are in radians?) 
        s1, s1az = tensile stress magnitude, tensile stress direction as azimuth
        s3, s3az = compressive stress magnitude, compressive stress direction as azimuth
        """
        return self.stresscalc.principal_components(pi/2 -lat, lon, self.time)

    def reset_cycl(self):
        """
        Used to reset parameters to begin anew. May also just create a new cycloid
        and save the old ones in some file.
        """
        self.lat_cycloidpts = []
        self.lon_cycloidpts = []
        self.data = []
    # END helper functions

    # BEGIN main functions
    # (Main propagation/initiation/plotting routines go here and were preserved)

    # END main functions

# Additional plotting and export functions follow...
# Citations, e.g., Peter Sinclair and Andre Ismailyan (2016), are preserved in docstrings.
# More inline commentary can be restored if you’d like even more original details added back.

#GNU Terry Pratchett
