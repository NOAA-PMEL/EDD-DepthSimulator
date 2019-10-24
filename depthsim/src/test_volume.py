import pytest
import numpy as np

from profiler import Volume


def test_class_volume_should_accept_direct_init():
    v = Volume(volume=32)

    assert 32 == v._volume

def test_class_volume_should_fail_if_volume_is_non_numeric():
    with pytest.raises(ValueError):
        v = Volume(volume='a')

def test_class_volume_should_get_volume_through_property():
    v = Volume(volume=13)

    assert 13 == v.volume

def test_class_volume_adding_together_creates_new_volume():
    v0 = Volume(volume=10)
    v1 = Volume(volume=20)
    v2 = v0 + v1

    assert 30 == v2.volume

def test_class_volume_subtracts_tegether_creates_new_volume():
    v0 = Volume(volume=100)
    v1 = Volume(volume=20)
    v2 = v0 - v1

    assert 80 == v2.volume

def test_class_volume_cannot_have_negative_volume():
    with pytest.raises(ValueError):
        v = Volume(volume=-50)