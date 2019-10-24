import pytest
import numpy as np

from profiler import Drag



# def test_drag_init_should_accept_shape_input():

#     d = Drag(shape='sphere')

#     assert 'sphere' == d._shape

def test_drag_init_should_fail_if_shape_input_invalid():
    with pytest.raises(ValueError):
        d = Drag(shape='Rabbit')

def test_drag_init_should_fail_if_shape_not_included():

    with pytest.raises(KeyError):
        d = Drag()

def test_drag_init_should_fail_if_circle_shape_doesnt_include_diameter():
    with pytest.raises(KeyError):
        d = Drag(shape='sphere')

def test_drag_init_should_accept_diameter_if_valid_shape():
    d = Drag(shape='sphere', diameter=5.0)

    area = np.pi * ((5.0/2) ** 2)
    assert d._area == area

def test_drag_init_should_fail_if_diameter_is_not_number():
    with pytest.raises(ValueError):
        d = Drag(shape='cone',diameter='a')

def test_drag_init_should_accept_length_of_side_for_non_spherical():
    d = Drag(shape='cube', length=4)

    area = 4**2

    assert d._area == area

def test_drag_init_should_fail_if_cube_length_is_not_number():
    with pytest.raises(ValueError):
        d = Drag(shape='cube', length='a')

def test_drag_init_should_fail_if_cube_doesnt_include_length():
    with pytest.raises(KeyError):
        d = Drag(shape='cube')

def test_drag_get_area_by_property():
    d = Drag(shape='sphere', diameter=4.2)

    expected_area = np.pi * ((4.2/2) ** 2)

    assert expected_area == d.area

def test_drag_get_drag_coefficient_by_property():
    d = Drag(shape='cube', length=14)

    assert 1.05 == d.drag_coefficient

def test_calculate_drag_should_calculate_correctly():
    density = 2
    velocity = 1
    coef = 0.5
    area = 10

    expected_drag = 5
    d = Drag(shape='cube', length=2)

    drag = d.calculate_drag(coef, density, area, velocity)

    assert expected_drag == drag 

def test_drag_should_report_zero_for_no_velocity():
    d = Drag(shape='sphere', diameter=10)

    assert 0.0 == d.drag

def test_drag_should_report_correctly_for_velocity():
    d = Drag(velocity=1.0, shape='sphere', diameter=5.0)

    expected_drag = 4729.5

    assert np.isclose(expected_drag, d.drag, rtol=1e0, atol=1e-02,equal_nan=False)


def test_drag_density_setter():
    d = Drag(shape='sphere', diameter=10)

    d.density = 2222.2

    assert 2222.2 == d.density

def test_drag_density_setter_raises_error_for_non_number():
    d = Drag(shape='sphere', diameter=10)

    with pytest.raises(ValueError):
        d.density = 'whoa'

def test_drag_velocity_setter():
    d = Drag(shape='sphere', diameter=10)

    d.velocity = 1144.22

    assert 1144.22 == d.velocity

def test_drag_velocity_setter_raises_error_for_non_number():
    d = Drag(shape='sphere', diameter=10)

    with pytest.raises(ValueError):
        d.velocity = 'now'