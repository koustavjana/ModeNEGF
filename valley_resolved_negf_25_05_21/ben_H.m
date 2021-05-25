N = 2*Ng+1*Nc;
Np = N;
O = zeros(N,N);
o = zeros(Np,Np);
H11 = zeros(N,N);
H22 = zeros(N,N);
H12 = t*eye(N);
H21 = t*eye(N);
h11 = zeros(Np,Np);
h22 = zeros(Np,Np);
h12 = t*eye(Np);
h21 = t*eye(Np);

Z11 = ones(1,Np);
Z22 = ones(1,Np);
for ii = 1:Np
    if mod(ii,2) == 0
        Z11(ii) = -1;
    end
    if mod(ii,2) ~= 0
        Z22(ii) = -1;
    end
end

Ugate = Vgg*[kind*ones(1,Ng) kind*linspace(1,-1,Nc) -kind*ones(1,Ng)];
Z11 = diag(Z11.*Ugate);
Z22 = diag(Z22.*Ugate);

for ii = 1:N-1
    if mod(ii,2) == 0 
        H11(ii,ii+1) = t;
        H11(ii+1,ii) = t;
        H22(ii,ii-1) = t;
        H22(ii-1,ii) = t;
    end
end
for ii = 1:Np-1
    if mod(ii,2) == 0 
        h11(ii,ii+1) = t;
        h11(ii+1,ii) = t;
        h22(ii,ii-1) = t;
        h22(ii-1,ii) = t;
    end
end

if mod(N,2) == 0
    H22(N,N-1) = t;
    H22(N-1,N) = t;
end
if mod(Np,2) == 0
    h22(Np,Np-1) = t;
    h22(Np-1,Np) = t;
end


salpha = [h11+Z11 h12;h21 h22+Z22];
sbeta_dag = [o h12; o o];
sbeta = [o o; h21 o];

alpha = [H11 H12;H21 H22];
beta_dag = [O H12; O O];
beta = [O O; H21 O];

sHk = @(kx) salpha + exp(1i*kx)*sbeta + exp(-1i*kx)*sbeta_dag;
Hk = @(kx) alpha + exp(1i*kx)*beta + exp(-1i*kx)*beta_dag;
vk = @(kx) 1i*exp(1i*kx)*beta - 1i*exp(-1i*kx)*beta_dag;
