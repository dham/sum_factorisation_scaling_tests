import tsfc
import coffee
import json
import numpy as np


def create_form(form_string, family, degree, dimension, operation):
    from firedrake import dx, inner, grad
    import firedrake as f

    f.parameters['coffee']['optlevel'] = 'O3'
    f.parameters['pyop2_options']['opt_level'] = 'O3'
    f.parameters['pyop2_options']['simd_isa'] = 'avx'
    f.parameters["pyop2_options"]["lazy_evaluation"] = False

    m = f.UnitSquareMesh(2, 2, quadrilateral=True)
    if dimension == 3:
        m = f.ExtrudedMesh(m, 2)

    fs = f.FunctionSpace(m, family, degree)
    v = f.TestFunction(fs)
    if operation == "matrix":
        u = f.TrialFunction(fs)
    elif operation == "action":
        u = f.Function(fs)
    else:
        raise ValueError("Unknown operation: %s" % operation)

    if form_string == "u * v * dx":
        return u * v * dx
    elif form_string == "inner(grad(u), grad(v)) * dx":
        return inner(grad(u), grad(v)) * dx

data = {}

for dimension in 2, 3:
    data[dimension] = {}

    for operation in ("action",): # ("matrix", "action"):
        data[dimension][operation] = {}
        for form_string in ("u * v * dx", "inner(grad(u), grad(v)) * dx"):
            print "%s %s in %d dimensions" % (form_string, operation, dimension)
            flops = {}
            for d in range(1, 10):
                form = create_form(form_string, "Q", d, dimension, operation)
                ast = tsfc.compile_form(form)[0].ast
                flops[d] = coffee.visitors.EstimateFlops().visit(ast)

                print d, flops[d]
                if d > 1:
                    print np.log(float(flops[d])/flops[d-1])/np.log(float(d)/float(d-1))
            data[dimension][operation][form_string] = flops

json.dump(data, file("data_count.json", "w"))
