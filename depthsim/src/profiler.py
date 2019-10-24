import math
import numpy as np
from scipy.constants import g

class Profiler:
    def __init__(self):
        pass



class Depth:
    def __init__(self, depth=0.0, density=1025):
        self._depth = depth
        self._pressure = 0.0
        self._density = density

    @property
    def depth(self):
        return self._depth

    @property
    def pressure(self):
        
        return self.calculate_pressure(self._depth, self._density)

    @staticmethod
    def calculate_pressure(depth, density):
        '''
        Calculate the static pressure of water column by depth and water density.  

        P = r * g * h

        (adapted from https://www.grc.nasa.gov/WWW/k-12/WindTunnel/Activities/fluid_pressure.html)

        where:
        P: is pressure (Pa)
        r(rho): is density of fluid
        g: acceleration of gravity (constant)
        h: height of fluid above object

        Args:
            depth: Depth of object (in meters)
            density: Density of water column (kg/m^3)
        Returns:
            pressure: Pressure on object (kPa)
        
        Raises:
            ValueError: Non-numeric entry for depth of density



        '''
        if not isinstance(depth, (int, float)):
            raise ValueError("Depth must be a numeric value")
        if not isinstance(density, (int, float)):
            raise ValueError("Density must be a numeric value")

        return density * depth * g / 1000


    @staticmethod
    def calculate_depth(pressure, density):
        '''
        Calculate the depth of object based on density and pressure.

        h = P / (r * g) 

        (adapted from https://www.grc.nasa.gov/WWW/k-12/WindTunnel/Activities/fluid_pressure.html)

        where:
        P: is pressure (Pa)
        r(rho): is density of fluid
        g: acceleration of gravity (constant)
        h: height of fluid above object

        Args:
            pressure: Presure on object (in kPa)
            density: Density of water column (kg/m^3)
        Returns:
            pressure: Depth of object (m)
        
        Raises:
            ValueError: Non-numeric entry for pressure of density

        '''

        if not isinstance(pressure, (int, float)):
            raise ValueError("Pressure must be a numeric value")
        if not isinstance(density, (int, float)):
            raise ValueError("Density must be a numeric value")

        return (pressure * 1000.0) / (density * g) 
        
class Buoyancy:
    def __init__(self, volume=1.0, density=1025.0):
        self._volume = volume
        self._density = density


    @property
    def density(self):
        return self._density

    @density.setter
    def density(self, value):
        if not isinstance(value, (int, float)):
            raise ValueError("Density must be a number")
        
        self._density = value

    @property
    def volume(self):
        return self._volume

    @property
    def buoyancy(self):
        self._buoyancy = self.calculate_buoyancy(self._volume, self._density)
        return self._buoyancy
    
    @staticmethod
    def calculate_buoyancy(volume, density):
        '''
        Calculate the bouyant force on the object based on volume and density

        Fb = V * r * g

        where:
        Fb: Buoyant force (kN)
        V: volume of object (m^3)
        r(rho): is desity of fluid (kg/m^3)
        g: acceleration of gravity (contant ~9.81m/s^2)

        Args:
            volume: Volume of object (m^3)
            density: Density of fluid (kg/m^3)
            
        Returns:
            Force: buoyant force (kN)
        
        Raises:
            ValueError: Non-numeric entry for volume of density
        '''
        if not isinstance(volume, (int, float)):
            raise ValueError("Volume must be a number")
        if not isinstance(density, (int, float)):
            raise ValueError("Density must be a number")
        return (density * volume * g) / 1000