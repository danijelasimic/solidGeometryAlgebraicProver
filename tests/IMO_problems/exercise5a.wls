poly = { x12*x13 + -1,
-x14,
-x15,
1*x8 + -x16,
x8*x12 + -x17,
-x7*x12 + -x18,
x7*x11 + -x8*x10 + -x19,
x8*x12 + -x20,
-x7*x12 + 1*x12 + -x21,
x7*x11 + -1*x11 + -x8*x10 + x8*1 + -x22,
-x23,
-1*x12 + -x24,
1*x11 + -x25,
x26*x15 + -x10*x15 + -x27*x14 + x11*x14,
x26*x16 + -x10*x16 + -x28*x14 + x12*x14,
x27*x16 + -x11*x16 + -x28*x15 + x12*x15,
x26*x18 + -1*x18 + -x27*x17,
x26*x19 + -1*x19 + -x28*x17,
x27*x19 + -x28*x18,
x29*x21 + -x30*x20,
x29*x22 + -x31*x20,
x30*x22 + -x31*x21,
x29*x24 + -x7*x24 + -x30*x23 + x8*x23,
x29*x25 + -x7*x25 + -x31*x23,
x30*x25 + -x8*x25 + -x31*x24,
x26 + -x29,
x27 + -x30,
x28 + -x31,
1^2 + -x32,
x7^2 + -2*x7*x10 + x10^2 + x8^2 + -2*x8*x11 + x11^2 + x12^2 + -x33,
1^2 + -2*1*x7 + x7^2 + x8^2 + -x34,
x10^2 + x11^2 + x12^2 + -x35,
x36 + -x32 + -x33,
x37 + -x34 + -x35}

variables = {x7, x8, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, zz, zz1}


condition = x34
poly = Insert[poly, condition* zz1 - 1,1]

statement1 = x36 + -x37
Print[GroebnerBasis[Insert[poly, statement1* zz - 1,1], variables]]

