import os
import pytest


@pytest.fixture(scope='session', autouse=True)
def root_directory():
    path_of_this_file = os.path.dirname(__file__)

    return path_of_this_file


@pytest.fixture(scope='session', autouse=True)
def files_directory(root_directory):
    return os.path.join(root_directory, 'files')
