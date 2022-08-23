function A = estimatePhase5pnt2(x2)

S1 = x2(1);
S2 = x2(2);
S3 = x2(3);
S4 = x2(4);
S5 = x2(5);
CR = 2*S3 - S5 - S1;
SR = 2*(S2 - S4);
A = atan2(CR,SR) + pi/2;