poly = { -x13,
-x14,
1*x8 + -x15,
x16,
x13*x17 + x14*x18 + x15*x19 + x16,
x17*x14 + -x10*x14 + -x18*x13 + x11*x13,
x17*x15 + -x10*x15 + -x19*x13 + x12*x13,
x18*x15 + -x11*x15 + -x19*x14 + x12*x14,
2*x20 + -x10 + -x17,
2*x21 + -x11 + -x18,
2*x22 + -x12 + -x19}

variables = {x7, x8, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, zz}

condition = -x15 + x8
Print[GroebnerBasis[Insert[poly, condition* zz - 1,1], variables]]

