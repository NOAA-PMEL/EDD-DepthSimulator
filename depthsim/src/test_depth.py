import pytest
import numpy as np

# from depth import Depth
from profiler import Depth

def test_depth_class_should_initialize_to_depth_zero():
    d = Depth()

    assert 0 == d._depth

def test_depth_class_should_read_depth_as_param():
    d = Depth()
    assert 0 == d.depth

def test_depth_class_calculate_pressure_should_fail_for_invalid_depth():
    d = Depth()

    with pytest.raises(ValueError):
        d.calculate_pressure('A', 20)

def test_depth_class_calculate_pressure_should_fail_for_invalid_density():
    d = Depth()

    with pytest.raises(ValueError):
        d.calculate_pressure(20, 'B')

def test_depth_class_should_calculate_pressure_from_depth():
    d = Depth()

    depth = 100.00 # m
    density = 1025.0 #kg/m^3
    expected_pressure = 1004.54 #kPa

    pressure = d.calculate_pressure(depth, density)

    assert np.isclose(expected_pressure, pressure, rtol=1e-03, atol=1e-05, equal_nan=False)


def test_depth_class_calculate_depth_should_fail_for_invalid_pressure():
    d = Depth()

    with pytest.raises(ValueError):
        d.calculate_depth('A',40)

def test_depth_class_calculate_depth_should_fail_for_invalid_density():
    d = Depth()

    with pytest.raises(ValueError):
        d.calculate_depth(40, 'C')

def test_depth_class_should_calculate_depth_from_pressure():
    d = Depth()

    pressure = 12000.0 # kPa
    density = 1025.0 # kg/m^3
    expected_depth = 1193.81 # m

    depth = d.calculate_depth(pressure, density)
    print(depth)
    assert np.isclose(expected_depth, depth, rtol=1e-03, atol=1e-05, equal_nan=False)


def test_depth_class_property_pressure():

    d = Depth(depth=250.0)

    expected_pressure = 2511.4
    pressure = d.pressure

    assert np.isclose(expected_pressure, pressure, rtol=1e-03, atol=1e-05, equal_nan=False)