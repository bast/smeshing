import os
import sys
from subprocess import Popen, PIPE
from cffi import FFI

ffi = FFI()

interface = Popen(['cc', '-E', os.path.join('distmesh', 'src', 'distance_functions.h')],
                  stdout=PIPE).communicate()[0].decode("utf-8")
ffi.cdef(interface)

if sys.platform == "darwin":
    suffix = 'dylib'
else:
    suffix = 'so'

lib = ffi.dlopen(os.path.join('build', 'libdistance_functions.{0}'.format(suffix)))
