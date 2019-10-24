import pytest
import numpy as np

from profiler import Buoyancy


def test_class_buoyancy_should_initialize_with_given_parameters():
    b = Buoyancy(volume=1.2, density=1024)
    assert 1024 == b._density
    assert 1.2 == b._volume

def test_class_buoyancy_should_return_density_through_property():
    b = Buoyancy()

    assert 1025.0 == b.density

def test_class_buoyancy_density_setter_should_fail_for_non_numeric():
    b = Buoyancy()

    with pytest.raises(ValueError):
        b.density = 'a'

def test_class_buoyancy_should_return_volume_through_property():
    b = Buoyancy()

    assert 1.0 == b.volume

def test_class_buoyancy_should_set_density_for_valid_number():
    b = Buoyancy()

    b.density = 1100.0

    assert 1100.0 == b.density

def test_calculate_buoyancy_should_fail_for_invalid_volume():
    b = Buoyancy()

    with pytest.raises(ValueError):
        b.calculate_buoyancy('a', 10)

def test_calculate_buoyancy_should_fail_for_invalid_density():
    b = Buoyancy()

    with pytest.raises(ValueError):
        b.calculate_buoyancy(100,'b')

def test_class_buoyancy_should_report_correct_buoyancy():
    b = Buoyancy()

    expected_buoyancy = 10052
    buoyancy = b.calculate_buoyancy(1.0, 1025)

    assert np.isclose(expected_buoyancy, buoyancy, rtol=1e-1, atol=1e-2, equal_nan=False)

def test_class_buoyancy_should_return_buoyancy_for_property():
    b = Buoyancy()

    expected_buoyancy = 10052

    assert np.isclose(expected_buoyancy, b.buoyancy, rtol=1e-03, atol=1e-05, equal_nan=False)
