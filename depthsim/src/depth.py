import math
import numpy as np
from scipy.constants import g

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