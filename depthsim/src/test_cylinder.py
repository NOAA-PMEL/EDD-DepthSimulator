import pytest
import numpy as np

from profiler import Cylinder

def test_cyclinder_class_should_accept_diameter_and_length():
    c = Cylinder(diameter=10, length=2)

    assert 10 == c.diameter
    assert 2 == c.length


def test_cylinder_class_should_fail_if_diameter_in_not_number():
    with pytest.raises(ValueError):
        c = Cylinder(diameter='a')

def test_cylinder_class_should_fail_if_length_in_not_number():
    with pytest.raises(ValueError):
        c = Cylinder(length='a')

def test_cylinder_class_should_calculate_area_correctly():
    c = Cylinder(diameter=10,length=2)  

    expected_area = 78.5398


    assert np.isclose(expected_area, c.area, rtol=1e-4, atol=1e-6, equal_nan= False)

def test_cylinder_class_should_calculate_volume_correctly():
    c = Cylinder(diameter=10, length=2)
    
    expected_volume = 157.0796

    assert np.isclose(expected_volume, c.volume, rtol=1e-4, atol=1e-6,equal_nan=False)

def test_two_cylinder_class_should_add_together():
    c0 = Cylinder(diameter=10, length=4)
    c1 = Cylinder(diameter=10, length=2)
    c = c0 - c1
    expected_volume = c0.volume - c1.volume

    assert expected_volume == c.volume

def test_calculate_cylinder_volume_works():
    c = Cylinder(diameter=10, length=2)
    expected_volume = 157.0796

    assert np.isclose(expected_volume, c.calculate_volume(10,2), rtol=1e-4, atol=1e-6,equal_nan=False)

def test_calculate_cylinder_area_works():
    
    diameter = 3

    expected_area = np.pi * (diameter/2)**2

    assert expected_area == Cylinder.calculate_area(diameter)