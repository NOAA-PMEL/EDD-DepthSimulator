import pytest
import numpy as np
from scipy.constants import g

from profiler import Profiler
from profiler import Cylinder
from profiler import Volume
from profiler import terminal_velocity



def test_terminal_velocity_calc_should_not_accept_invalid_density():
    with pytest.raises(ValueError):
        terminal_velocity('density', 1, 2, 3)

def test_terminal_velocity_calc_should_not_accept_invalid_mass():
    with pytest.raises(ValueError):
        terminal_velocity(1, 'mass', 2, 3)

def test_terminal_velocity_calc_should_not_accept_invalid_coef():
    with pytest.raises(ValueError):
        terminal_velocity(1, 2, 'coef', 3)

def test_terminal_velocity_calc_should_not_accept_invalid_area():
    with pytest.raises(ValueError):
        terminal_velocity(1, 2, 3, 'area')

def test_terminal_velocity_calc_should_return_valid_answer():
    density = 1.5
    mass = 10
    area = 32
    coef = 0.82
    expected_vt = 2.23228
    vt = terminal_velocity(density, mass, coef, area)

    assert np.isclose(expected_vt, vt, rtol = 1e-3, atol=1e-4, equal_nan=False) 


def test_class_profiler_should_initialize_with_given_parameters():
    parameters = {'body_d': 0.4,
                    'body_l':1.0,
                    'piston_d':0.1,
                    'piston_l':0.3,
                    'density':1023.2,
                    'depth':0.0,
                    'velocity':0.0,
                    'mass':12.2
    }

    p = Profiler(**parameters)

    expected_body_volume = 0.12566
    expected_piston_volume = 0.002356

    assert np.isclose(expected_body_volume, p.body.volume, rtol=1e-3,atol=1e-4)
    assert np.isclose(expected_piston_volume, p.piston.volume, rtol=1e-3,atol=1e-4)
    assert parameters['density'] == p.water.density
    assert parameters['depth'] == p.water.depth
    assert parameters['velocity'] == p.velocity
    assert parameters['mass'] == p._mass

def test_class_profiler_should_calculate_total_volume_correctly():
    parameters = {'body_d': 0.4,
                    'body_l':1.0,
                    'piston_d':0.1,
                    'piston_l':0.3,
                    'density':1023.2,
                    'depth':0.0,
                    'velocity':0.0,
                    'mass':12.2
    }

    p = Profiler(**parameters)

    expected_body_volume = 0.12566
    expected_piston_volume = 0.002356
    expected_total_volume = expected_body_volume + expected_piston_volume

    assert np.isclose(expected_total_volume, p.volume, rtol=1e-4, atol=1e-6)


def test_class_profiler_should_calculate_total_profiler_density_correctly():
    parameters = {'body_d': 0.4,
                    'body_l':1.0,
                    'piston_d':0.1,
                    'piston_l':0.3,
                    'density':1023.2,
                    'depth':0.0,
                    'velocity':0.0,
                    'mass':12.2
    }

    p = Profiler(**parameters)

    expected_body_volume = 0.12566
    expected_piston_volume = 0.002356
    expected_total_volume = expected_body_volume + expected_piston_volume
    expected_profiler_density = parameters['mass'] / expected_total_volume

    assert np.isclose(expected_profiler_density, p.density, rtol=1e-4, atol=1e-6,equal_nan=False)


def test_class_profiler_should_calculate_total_acceleration_correctly_balanced_forces():
    # Testing to see acceleration = 0 when F_gravity == F_buoyancy

    parameters = {
        'body_d': 0.50465,
        'piston_d':0.0,
        'piston_l':0.0,
        'density':1025,
        'depth':0.0,
        'velocity':0.0,
        'mass':50.0
    }

    # f_b  - f_g = 0
    # f_g = m * g
    # f_b = V * r * g

    # f_b = m * acc_b
    # m * acc_b - m * g = 0
    # acc_b - g = 0
    f_g = parameters['mass'] * g

    volume = parameters['mass'] / parameters['density']
    area = ((parameters['body_d']/2)**2) * np.pi
    length = volume/area

    parameters['body_l'] = length

    p = Profiler(**parameters)

    assert 0 == p.acceleration

def test_class_profiler_should_calculate_total_acceleration_with_drag_included():
    parameters = {
        'body_d': 0.50465,
        'body_l': 1.0,
        'piston_d':0.0,
        'piston_l':0.0,
        'density':1025,
        'depth':0.0,
        'velocity':1.0,
        'mass':50.0
    }

    expected_f_drag = 84.05

    f_buoy = 2010.55
    f_grav = 490.33

    f_total = f_buoy - (expected_f_drag + f_grav)

    p = Profiler(**parameters)
    print(p.acceleration)
    assert np.isclose(expected_f_drag, p.drag.drag, rtol=1e-4, atol=1e-6, equal_nan=False)