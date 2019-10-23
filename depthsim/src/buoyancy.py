import math
import numpy as np
from scipy.constants import g

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


