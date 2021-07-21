import os
import pytest
from cantherm import get_sample_file_path


@pytest.mark.parametrize("filename", [("bz.log"), ("phosphonyl.log")])
def test_get_data(filename):
    path = get_sample_file_path(filename)
    assert os.path.exists(path)


@pytest.mark.parametrize("filename", [("fake1.loggg"), ("fake2.ou")])
def test_get_data_fail(filename):
    with pytest.raises(ValueError):
        get_sample_file_path(filename)
