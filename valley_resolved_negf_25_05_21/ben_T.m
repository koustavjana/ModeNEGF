%%
dim = 2*NL*Np+4*N;
H = zeros(dim,dim);
Uself = pot*eye(dim);
Uself(1:2*N,1:2*N) = zeros(2*N,2*N);
Uself(dim-2*N+1:dim,dim-2*N+1:dim) = zeros(2*N,2*N);
H(1:2*N,1:2*N) = alpha;
H(dim-2*N+1:dim,dim-2*N+1:dim) = alpha;
H(N+1:N+Np,1+2*N:2*N+Np) = h21;
H(1+2*N:2*N+Np,N+1:N+Np) = h12;
H(dim-2*N+1:dim-2*N+Np,dim-2*N-Np+1:dim-2*N) = h12;
H(dim-2*N-Np+1:dim-2*N,dim-2*N+1:dim-2*N+Np) = h21;

for ii = 1:NL
    H(1+2*N+2*(ii-1)*Np:2*N+2*ii*Np,1+2*N+2*(ii-1)*Np:2*N+2*ii*Np) = salpha;
end
for ii = 1:NL-1
    H(1+2*N+2*(ii-1)*Np:2*N+2*ii*Np,1+2*N+2*(ii)*Np:2*N+2*(ii+1)*Np) = sbeta;
end
for ii = 2:NL
    H(1+2*N+2*(ii-1)*Np:2*N+2*ii*Np,1+2*N+2*(ii-2)*Np:2*N+2*(ii-1)*Np) = sbeta_dag;
end

HH = H;
HH(1:N,1:N) = H(1:N,1:N)+ML;
HH(1+dim-N:dim,1+dim-N:dim) = H(1+dim-N:dim,1+dim-N:dim)+MR;

GR = inv(E*eye(dim) - HH - Uself);
GA = GR';

TL = [TL O;O O];
TR = [O O; O TR];
ML = [ML O;O O];
MR = [O O; O MR];

ML_dum = zeros(dim,dim);
MR_dum = zeros(dim,dim);
ML_dum(1:2*N,1:2*N) = ML;
MR_dum(dim-2*N+1:dim,dim-2*N+1:dim) = MR;
ML = ML_dum;
MR = MR_dum;
TL = 1i*(ML-ML');
TR = 1i*(MR-MR');

% Trans = trace(TL*GR*TR*GA)
% Trans = trace(TR*GR*TL*GA)