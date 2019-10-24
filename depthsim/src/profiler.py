import math
import numpy as np
from scipy.constants import g


class Profiler:
    def __init__(self, **kwargs):
        self.body = Cylinder(kwargs["body_d"], kwargs["body_l"])
        self.piston = Cylinder(kwargs["piston_d"], kwargs["piston_l"])
        self.water = Water(kwargs["density"], kwargs["depth"])
        self.depth = Depth(self.water.depth, self.water.density)
        self.buoyancy = Buoyancy(self.volume, self.water.density)
        self.velocity = kwargs["velocity"]
        self.drag = Drag(
            self.water.density,
            self.velocity,
            shape="long_cylinder",
            diameter=self.body.diameter,
        )
        self._mass = kwargs["mass"]

    @property
    def volume(self):
        return self.body.volume + self.piston.volume

    @property
    def velocity(self):
        return self._velocity

    @velocity.setter
    def velocity(self, value):
        self._velocity = value

    @property
    def density(self):
        return self._mass / self.volume

    @property
    def acceleration(self):
        self.buoyancy.volume = self.volume
        f_buoy = self.buoyancy.buoyancy
        print(f'f_buoy= {f_buoy}')
        f_grav = self._mass * g
        print(f'f_grav= {f_grav}')
        f_drag = self.drag.drag 
        print(f'f_drag = {f_drag}')
        accel = 0
        if f_buoy > f_grav:
            accel = f_buoy - (f_grav + f_drag)
        if f_buoy < f_grav:
            accel = (f_buoy + f_drag) - f_grav

        return accel

class Water:
    def __init__(self, density=1025, depth=0.0):
        self._density = density
        self._depth = depth

    @property
    def density(self):
        return self._density

    @density.setter
    def density(self, value):
        self._density = value

    @property
    def depth(self):
        return self._depth

    @depth.setter
    def depth(self, value):
        self._depth = value


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
        """
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



        """
        if not isinstance(depth, (int, float)):
            raise ValueError("Depth must be a numeric value")
        if not isinstance(density, (int, float)):
            raise ValueError("Density must be a numeric value")

        return density * depth * g / 1000

    @staticmethod
    def calculate_depth(pressure, density):
        """
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

        """

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

    @volume.setter
    def volume(self, value):
        self._volume = value

    @property
    def buoyancy(self):
        self._buoyancy = self.calculate_buoyancy(self._volume, self._density)
        return self._buoyancy

    @staticmethod
    def calculate_buoyancy(volume, density):
        """
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
        """
        if not isinstance(volume, (int, float)):
            raise ValueError("Volume must be a number")
        if not isinstance(density, (int, float)):
            raise ValueError("Density must be a number")
        return (density * volume * g)


SHAPE_DRAG_COEF = {
    "sphere": 0.47,
    "hemisphere": 0.42,
    "cone": 0.50,
    "cube": 1.05,
    "angled_cube": 0.80,
    "long_cylinder": 0.82,
    "short_cylinder": 1.15,
    "streamlined_body": 0.04,
    "streamlined_half": 0.09,
}


class Drag:
    def __init__(self, density=1025, velocity=0.0, **kwargs):

        self.velocity = velocity
        self.density = density

        if "shape" not in kwargs.keys():
            raise KeyError("Missing Shape of Profiler")

        if kwargs["shape"].lower() in SHAPE_DRAG_COEF:
            self._shape = kwargs["shape"].lower()
        else:
            raise ValueError(f"Invalid shape: {kwargs['shape']}")

        if (
            (self._shape == "sphere")
            or (self._shape == "hemisphere")
            or (self._shape == "cone")
            or (self._shape == "long_cylinder")
            or (self._shape == "short_cylinder")
        ):

            if "diameter" not in kwargs.keys():
                raise KeyError("Diameter required")
            if not isinstance(kwargs["diameter"], (int, float)):
                raise ValueError("Diameter must be a number")

            self._area = np.pi * ((kwargs["diameter"] / 2) ** 2)

        elif self._shape == "cube":
            if "length" not in kwargs.keys():
                raise KeyError("Length required")
            if not isinstance(kwargs["length"], (int, float)):
                raise ValueError("Length must be a number")
            self._area = kwargs["length"] ** 2
        else:
            raise TypeError("What is it?")

        self._drag_coefficient = SHAPE_DRAG_COEF[self._shape]

        self.drag

    @property
    def area(self):
        return self._area

    @property
    def drag_coefficient(self):
        return self._drag_coefficient

    @property
    def density(self):
        return self._density

    @property
    def velocity(self):
        return self._velocity

    @velocity.setter
    def velocity(self, value):
        if not isinstance(value, (int, float)):
            raise ValueError("Velocity must be a number")
        self._velocity = value

    @density.setter
    def density(self, value):
        if not isinstance(value, (int, float)):
            raise ValueError("Density must be a number")
        self._density = value

    @property
    def drag(self):
        return self.calculate_drag(
            self.drag_coefficient, self.density, self.area, self.velocity
        )

        return self._drag

    @staticmethod
    def calculate_drag(drag_coef, density, area, velocity):
        """ 
        Calculate the drag force (N) of an object

        Fd = 1/2 * r * u^2 * Cd * A

        where:

        Fd: Drag force (N)
        r(rho): mass density of the fluid (kg/m^3)
        u: flow velocity relative to object (m/s)
        A: reference area (m^2)
        Cd: Drag coefficient (dimensionless)

        Args:
            drag_coef: Drag coefficent (Cd)
            density: Mass density (rho)
            area: Reference area (A)
            velocity: Flow velocity (u)
        
        Returns:
            Drag force: Fd
        """
        return (density * (velocity ** 2) * drag_coef * area) / 2


class Volume:
    def __init__(self, **kwargs):
        if "volume" in kwargs.keys():
            if not isinstance(kwargs["volume"], (int, float)):
                raise ValueError("Volume must be a number")
            self.volume = kwargs["volume"]

    @property
    def volume(self):
        return self._volume

    @volume.setter
    def volume(self, value):
        if value < 0:
            raise ValueError("Volume must be a positive number")
        self._volume = value

    def __add__(self, other):
        return Volume(volume=(self.volume + other.volume))

    def __sub__(self, other):
        return Volume(volume=(self.volume - other.volume))


class Cylinder(Volume):
    def __init__(self, diameter=0, length=0):
        if not isinstance(diameter, (int, float)):
            raise ValueError("diameter must be a number")
        if not isinstance(length, (int, float)):
            raise ValueError("length must be a number")

        self._diameter = diameter
        self._length = length

    @property
    def diameter(self):
        return self._diameter

    @property
    def length(self):
        return self._length

    @property
    def area(self):
        return ((self._diameter / 2) ** 2) * np.pi

    @property
    def volume(self):
        return self.calculate_volume(self.diameter, self.length)
        # return self.area * self._length

    @staticmethod
    def calculate_area(diameter):
        return (diameter / 2) ** 2 * np.pi

    @staticmethod
    def calculate_volume(diameter, length):
        return ((diameter / 2) ** 2) * np.pi * length


def terminal_velocity(density_f, mass, drag_coef, area):
    """ 
    Calculates the terminal velocity of an object moving through fluid

    Vt = sqrt( (2 * m * g) / (r * A * Cd) )
    
    where:
        m: Mass of object 
        g: gravitational constant 
        r(rho): density of fluid
        A: projected area
        Cd: coefficient of drag

    Args:
        density_f: Fluid density (kg/m^3)
        mass: Mass of object (kg)
        drag_coef: Cofficient of drag for object
        area: Project area of the object

    Returns:
        vt: Terminal Velocity of object
    
    Raises:
        ValueError: For invalid entries

    """
    if not isinstance(density_f, (int, float)):
        raise ValueError("density_f must be a number")

    if not isinstance(mass, (int, float)):
        raise ValueError("mass must be a number")
    if not isinstance(drag_coef, (int, float)):
        raise ValueError("drag_coef must be a number")

    if not isinstance(area, (int, float)):
        raise ValueError("area must be a number")

    return np.sqrt((2 * mass * g) / (density_f * area * drag_coef))
