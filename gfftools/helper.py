import logging
from typing import Union, List, Tuple

log = logging.getLogger(__name__)


def _get(index: int, list: Union[List, Tuple], default):
    """
    Get an element from a list or return a default.
    """
    try:
        return list[index]
    except IndexError:
        return default


def _get_and_cast(index: int, list: Union[List, Tuple], type_def, default):
    """
    Get an element from a list and cast. Return a default if item does not exist.
    """

    try:
        element = list[index]
    except IndexError:
        return default

    try:
        return type_def(element)
    except TypeError:
        return default