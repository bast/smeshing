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


def get_resolution(points,
                   use_tanh,
                   ref_points,
                   nearest_distance_at_coastline_point,
                   flanders_indices):

    num_points = len(points)
    (x, y) = zip(*points)

    num_reference_points = len(ref_points)
    (points_x, points_y) = zip(*ref_points)

    resolutions_np = np.zeros(num_points, dtype=np.float64)
    resolutions_p = _ffi.cast("double *", resolutions_np.ctypes.data)

    _lib.get_resolution(num_points,
                        x,
                        y,
                        resolutions_p,
                        use_tanh,
                        num_reference_points,
                        points_x,
                        points_y,
                        nearest_distance_at_coastline_point,
                        flanders_indices)

    return resolutions_np.tolist()
