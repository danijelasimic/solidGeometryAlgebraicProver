poly = { 2*x25 + -1,
2*x26,
2*x27,
2*x28,
2*x29 + -1,
2*x30,
x31 + -x25 + x28,
x32 + -x26 + x29,
x33 + -x27 + x30,
-x34,
-x35,
1*1 + -x36,
x37,
x38*x34 + x39*x35 + x40*x36,
x38*x25 + x39*x26 + x40*x27 + x41,
x31*x38 + x32*x39 + x33*x40,
x43*1 + -x43*1 + -x44*1 + 1*1,
x44*1 + -1*1 + -x42*1 + x42*1,
x42*1 + -x43*1,
x38*x42 + x39*x43 + x40*x44 + x41,
x45 + -x42^2 + -x43^2 + -1^2 + 2*1*x44 + -x44^2,
x46 + -x25^2 + 2*x25*x28 + -x28^2 + -x26^2 + 2*x26*x29 + -x29^2 + -x27^2 + 2*x27*x30 + -x30^2}

variables = {x25, x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39, x40, x41, x42, x43, x44, x45, x46, zz, zz1}

statement1 = 4*x45 + -x46
poly = Insert[poly, statement1* zz - 1,1]
poly = Insert[poly, (-1 + 2*x46)*zz1 - 1, 1]
Print[GroebnerBasis[poly, variables]]

