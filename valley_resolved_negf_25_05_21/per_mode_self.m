clear;
clc;
% close all;

Ng = 12;
Nc = 4;
NL = 17;
t = -1; 
Vgg = 1;
kind = 1;
E_arr = 0.5;
pot_arr = E_arr-Vgg:Vgg/49:E_arr+Vgg;
% pot_arr = E_arr;

Kpval = [];
Knval = [];
for E = E_arr
for pot = pot_arr 
ben_H
% ben_Ek
ben_M
ben_self
ben_T

% H12*inv(E*eye(N)-H22)*H21*(inv(Fn1)+eye(N))
ML_mode = [];
MR_mode = [];

Un1_inv = inv(Un1);
Up2_inv = inv(Up2);
for ii = 1:N
    ML_mode = [ML_mode H12*inv(E*eye(N)-H22)*H21*((Un1(:,ii)*(1/Yn(ii,ii))*Un1_inv(ii,:))+(Un1(:,ii)*Un1_inv(ii,:)))];
    MR_mode = [MR_mode H21*inv(E*eye(N)-H11)*H12*((Up2(:,ii)*(Yp(ii,ii))*Up2_inv(ii,:))+(Up2(:,ii)*Up2_inv(ii,:)))];
end

for ii = 1:N 
    for jj = 1:N 
        TL_mode = 1i*(ML_mode(:,1+(jj-1)*N:jj*N)-(ML_mode(:,1+(jj-1)*N:jj*N))');
        TR_mode = 1i*(MR_mode(:,1+(ii-1)*N:ii*N)-(MR_mode(:,1+(ii-1)*N:ii*N))');
        T(ii,jj) = trace(TR_mode*GR(dim-N+1:dim,1:N)*TL_mode*GA(1:N,dim-N+1:dim));     
    end
end

Tkp = 0;
Tkn = 0;
for ii = 1:N
    TV(ii) = sum(T(ii,:),'all');
    if abs(abs(Yp(ii,ii))-1) < 1e-6 
    if real(log(Yp(ii,ii))/1i) > 0
        Tkp = Tkp+TV(ii);
    else
        Tkn = Tkn+TV(ii);
    end
    end
end

Kpval = [Kpval real(Tkp)];
Knval = [Knval real(Tkn)];
Tkp
Tkn

end
end

% plot((E_arr-pot_arr)./del,Kpval,'LineWidth',2);
% hold on
figure(1)
plot((E_arr-pot_arr),(Kpval-Knval)./(Kpval+Knval),'LineWidth',2);
hold off
