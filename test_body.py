import firedrake as f

f.parameters["pyop2_options"]["lazy_evaluation"] = False


def form_assembler(form, n, family, degree, dimension, operation):

    m = f.UnitSquareMesh(n, n, quadrilateral=True)
    if dimension == 3:
        m = f.ExtrudedMesh(m, n)

    fs = f.FunctionSpace(m, family, degree)
    v = f.TestFunction(fs)
    if operation == "matrix":
        u = f.TrialFunction(fs)
    elif operation == "action":
        u = f.Function(fs)
    else:
        raise ValueError("Unknown operation: %s" % operation)

    actualform = form(u, v)

    if operation == "matrix":
        return lambda: f.assemble(actualform).M
    else:
        return lambda: f.assemble(actualform)
