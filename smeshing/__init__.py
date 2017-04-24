import os
from cffi import FFI
from .cffi_helpers import get_lib_handle
import numpy as np


_this_path = os.path.dirname(os.path.realpath(__file__))

_build_dir = os.getenv('SMESHING_BUILD_DIR')
if _build_dir is None:
    _build_dir = _this_path
else:
    _build_dir = os.path.join(_build_dir, 'lib')

_include_dir = _this_path

_lib = get_lib_handle(
    ['-DSMESHING_API=', '-DCPP_INTERFACE_NOINCLUDE'],
    'smeshing.h',
    'smeshing',
    _build_dir,
    _include_dir
)

_ffi = FFI()


def get_resolution(x,
                   y,
                   use_tanh,
                   points,
                   nearest_distance_at_coastline_point,
                   flanders_indices):

    num_points = len(points)
    (points_x, points_y) = zip(*points)

    return _lib.get_resolution(x,
                               y,
                               use_tanh,
                               num_points,
                               points_x,
                               points_y,
                               nearest_distance_at_coastline_point,
                               flanders_indices)
