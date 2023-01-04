import os, numpy

numpy_include_dir = numpy.get_include()
include_file = os.path.join(numpy_include_dir, "numpy", "arrayobject.h").replace("\\", "/")
with open("numpy_arrayobject.h", "wt") as f:
    f.write(f"#include \"{include_file}\"")
print(numpy_include_dir, end='')