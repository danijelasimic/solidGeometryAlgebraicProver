poly = { x7 + -x10 + -x4 + x1,
x8 + -x11 + -x5 + x2,
x9 + -x12 + -x6 + x3,
x19,
x20 + -x17,
x21 + -x18 + x15,
x25*x27 + x26*x27 + -x26*x1,
x25*x28 + x26*x28 + -x26*x2,
x25*x29 + x26*x29 + -x26*x3 + -x25*x15,
x25*x30 + x26*x30 + -x26*x4,
x25*x31 + x26*x31 + -x26*x5 + -x25*x17,
x25*x32 + x26*x32 + -x26*x6 + -x25*x18,
x25*x33 + x26*x33 + -x26*x7 + -x25*x19,
x25*x34 + x26*x34 + -x26*x8 + -x25*x20,
x25*x35 + x26*x35 + -x26*x9 + -x25*x21,
x25*x36 + x26*x36 + -x26*x10,
x25*x37 + x26*x37 + -x26*x11,
x25*x38 + x26*x38 + -x26*x12}

variables = {x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x15, x17, x18, x19, x20, x21, x25, x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, zz, zz1}

condition = x25 + x26
poly = Insert[poly, condition* zz1 - 1,1]

statement1 = x31*x38 + -x31*x35 + -x28*x38 + x28*x35 + -x32*x37 + x32*x34 + x29*x37 + -x29*x34
Print[GroebnerBasis[Insert[poly, statement1* zz - 1,1], variables]]

statement2 = x32*x36 + -x32*x33 + -x29*x36 + x29*x33 + -x30*x38 + x30*x35 + x27*x38 + -x27*x35
Print[GroebnerBasis[Insert[poly, statement2* zz - 1,1], variables]]

statement3 = x30*x37 + -x30*x34 + -x27*x37 + x27*x34 + -x31*x36 + x31*x33 + x28*x36 + -x28*x33
Print[GroebnerBasis[Insert[poly, statement3* zz - 1,1], variables]]

