from firedrake import *
parameters['coffee']['optlevel'] = 'O3'
parameters['pyop2_options']['opt_level'] = 'O3'
parameters['pyop2_options']['simd_isa'] = 'avx'

m = UnitSquareMesh(2, 2, quadrilateral=True)

fs = FunctionSpace(m, "Q", 5)

u = TrialFunction(fs)
v = TestFunction(fs)

print tsfc_interface.compile_form(u * v * dx, "foo")
