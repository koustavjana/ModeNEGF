%% Topo valley
clear;
clc;
close all;

Ng = 12;
Nc = 4;
NL = 17;
t = -1; 
Vgg = 1;
kind = 1;
% del = -pi*t/40;
% E_arr = -0.88:0.0056:0.88;
E_arr = 0.5;
% pot_arr = E_arr-4*del:del/49:E_arr+4*del;
% pot_arr = 0;
pot_arr = E_arr-Vgg:Vgg/49:E_arr+Vgg;
% pot_arr = E_arr;
fL = 1;
fR = 1;

NEGF = [];
Scatter = [];
Kpval = [];
Knval = [];

for E = E_arr
for pot = pot_arr    
ben_H
% ben_Ek
ben_M
ben_self
ben_T

MLin = TL*fL;
MRin = TR*fR;
Min = TL*fL + TR*fR;
Gn = GR*Min*GA;

Uk = [Ukp Ukn];
Yk = diag(inv(Uk)*M*Uk);
Kk = log(Yk)/1i;
vel = [];
velp = [];
veln = [];
Trans = trace(TR*GR*TL*GA);
disp(Trans)
NEGF = [NEGF real(Trans)];

for jj = 1:2*N
    vop = vk(Kk(jj));
    vel(jj) = (Uk(:,jj)')*vop*Uk(:,jj);
end
for jj = 1:N
    vop = vk(log(Yp(jj,jj))/1i);
    velp(jj) = (Up(:,jj)')*vop*Up(:,jj);
    vop = vk(log(Yn(jj,jj))/1i);
    veln(jj) = ((Un(:,jj)')*vop*Un(:,jj));
end

Uk_inv = inv(Uk);
gL_inv = H12*inv(E*eye(N)-H22)*H21*(inv(Fp1)-inv(Fn1));
gR_inv = H21*inv(E*eye(N)-H11)*H12*(Fn2-Fp2);
TT = inv(Up2)*GR(dim-N+1:dim,1:N)*gL_inv*Up1;
TT_r = inv(Un1)*GR(1:N,dim-N+1:dim)*gR_inv*Un2;


VR_2 = zeros(size(TT));
VL_inv_2 = zeros(size(TT));
VL_2 = zeros(size(TT));
VR_inv_2 = zeros(size(TT));

for ii = 1:N     
if abs(abs(Yp(ii,ii))-1) < 1e-6    
   VR_2(ii,ii) = sqrt(abs(velp(ii)));
   VL_inv_2(ii,ii) = sqrt(1/abs(velp(ii)));
end
end

for ii = 1:N     
if abs(abs(Yn(ii,ii))-1) < 1e-6    
   VL_2(ii,ii) = sqrt(abs(veln(ii)));
   VR_inv_2(ii,ii) = sqrt(1/abs(veln(ii)));
end
end

tau = VR_2*TT*VL_inv_2;
tau_r = VL_2*TT_r*VR_inv_2;

T = tau*tau';
T_r = tau_r'*tau_r;
disp(trace(T))
disp(trace(T_r))

Scatter = [Scatter trace(T)];

Tkp = 0;
Tkn = 0;
Tkp_r = 0;
Tkn_r = 0;

for ii = 1:N
    TV(ii) = T(ii,ii);
    if real(log(Yp(ii,ii))/1i) > 0
        Tkp = Tkp+TV(ii);
    else
        Tkn = Tkn+TV(ii);
    end
end
for ii = 1:N
    TV_r(ii) = T_r(ii,ii);
    if real(log(Yn(ii,ii))/1i) > 0
        Tkp_r = Tkp_r+TV_r(ii);
    else
        Tkn_r = Tkn_r+TV_r(ii);
    end
end
Kpval = [Kpval Tkp];
Knval = [Knval Tkn];
Tkp
Tkn
% Tkp_r
% Tkn_r




end
end


plot((E_arr-pot_arr),(Kpval-Knval)./(Kpval+Knval),'LineWidth',2);
hold on



% figure(2)
% % plot(E_arr,NEGF);
% % hold on
% plot((E_arr-pot_arr)./del,Scatter,'LineWidth',2);
% hold on
% % plot((E_arr-pot_arr)./del,Kpval,'LineWidth',2);
% % hold on
% % plot((E_arr-pot_arr)./del,Knval,'LineWidth',2);
% % hold off
