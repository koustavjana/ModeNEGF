ML = H12*inv(E*eye(N)-H22)*H21*(inv(Fn1)+eye(N));
MR = H21*inv(E*eye(N)-H11)*H12*(Fp2+eye(N));
TL = 1i*(ML-ML');
TR = 1i*(MR-MR');