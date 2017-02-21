import timeit
import json
import numpy as np


def setup(form, n, family, degree, dimension, operation):

    return """
import test_body
from firedrake import dx, inner, grad

form = lambda u, v: %(form)s

a = test_body.form_assembler(form, %(n)d, "%(family)s", %(degree)d, %(dimension)d, "%(operation)s")
""" % {"form": form,
       "n": n,
       "family": family,
       "degree": degree,
       "dimension": dimension,
       "operation": operation}

data = {}

for dimension in 2, 3:
    data[dimension] = {}
    size = 10 ** (6 / dimension)  # Size along one dimension of the mesh. Produce a million cell mesh at degree 1.
    for operation in ("action",): #("matrix", "action"):
        data[dimension][operation] = {}
        for form in ("u * v * dx", "inner(grad(u), grad(v)) * dx"):
            print "%s %s in %d dimensions" % (form, operation, dimension)
            times = {}
            for d in range(1, 10):
                scaled_size = size // d  # Attempt to keep the memory footprint approximately constant.
                size_scale = float(d) / size
                raw_times = timeit.repeat("a()", setup(form, scaled_size, "Q", d, dimension, operation), number=1)
                times[d] = [t * size_scale ** dimension for t in raw_times]
                print d, times[d]
                if d > 1:
                    print np.log(min(times[d])/min(times[d-1]))/np.log(float(d)/float(d-1))
            data[dimension][operation][form] = times

json.dump(data, file("data.json", "w"))
