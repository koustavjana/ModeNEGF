clear;
clc;
% close all;

Ng = 2;
Nc = 2;
NL = 10;
t = -1; 
Vgg = 0.0;
kind = -1;
E_arr = 0.9;
pot_arr = E_arr-Vgg*1.5:Vgg/50:E_arr+Vgg*1.5;
pot_arr = 0;

fL = 1;
fR = 0;

Kpval = [];
Knval = [];
for E = E_arr
for pot = pot_arr 
ben_H
ben_Ek
ben_M
ben_self
ben_T

MLin = TL*fL;
MRin = TR*fR;
Min = TL*fL + TR*fR;
Gn = GR*Min*GA;
A = 1i*(GR-GA);
Gp = A - Gn;

ML_mode = [];
MR_mode = [];
MLin_mode = [];
MRin_mode = [];

Un1_inv = inv(Un1);
Up2_inv = inv(Up2);
Up1_inv = inv(Up1);
Un2_inv = inv(Un2);
for ii = 1:N
    ML_mode = [ML_mode H12*inv(E*eye(N)-H22)*H21*((Un1(:,ii)*(1/Yn(ii,ii))*Un1_inv(ii,:))+(Un1(:,ii)*Un1_inv(ii,:)))];
    MR_mode = [MR_mode H21*inv(E*eye(N)-H11)*H12*((Up2(:,ii)*(Yp(ii,ii))*Up2_inv(ii,:))+(Up2(:,ii)*Up2_inv(ii,:)))];
    MLin_mode = [MLin_mode H12*inv(E*eye(N)-H22)*H21*((Up1(:,ii)*(Yp(ii,ii))*Up1_inv(ii,:))+(Up1(:,ii)*Up1_inv(ii,:)))];
    MRin_mode = [MRin_mode H21*inv(E*eye(N)-H11)*H12*((Un2(:,ii)*(1/Yn(ii,ii))*Un2_inv(ii,:))+(Un2(:,ii)*Un2_inv(ii,:)))];
end

for ii = 1:N 
    for jj = 1:N 
        TL_mode = 1i*(MLin_mode(:,1+(jj-1)*N:jj*N)-(MLin_mode(:,1+(jj-1)*N:jj*N))');
        TR_mode = 1i*(MR_mode(:,1+(ii-1)*N:ii*N)-(MR_mode(:,1+(ii-1)*N:ii*N))');
        T(ii,jj) = trace(TR_mode*GR(dim-N+1:dim,1:N)*TL_mode*GA(1:N,dim-N+1:dim));     
    end
end

Tkp = 0;
Tkn = 0;
MoutL_kp = zeros(size(H11));
MoutL_kn = zeros(size(H11));
MoutR_kp = zeros(size(H11));
MoutR_kn = zeros(size(H11));
MinL_kp = zeros(size(H11));
MinL_kn = zeros(size(H11));
MinR_kp = zeros(size(H11));
MinR_kn = zeros(size(H11));


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

for ii = 1:N
    if abs(abs(Yp(ii,ii))-1) < 1e-6 
    if real(log(Yp(ii,ii))/1i) > 0
        MoutR_kp = MoutR_kp + (1-fR)*1i*(MR_mode(:,1+(ii-1)*N:ii*N)-(MR_mode(:,1+(ii-1)*N:ii*N))');
        MinL_kp = MinL_kp + fL*1i*(MLin_mode(:,1+(ii-1)*N:ii*N)-(MLin_mode(:,1+(ii-1)*N:ii*N))');
    else
        MoutR_kn = MoutR_kn + (1-fR)*1i*(MR_mode(:,1+(ii-1)*N:ii*N)-(MR_mode(:,1+(ii-1)*N:ii*N))');
        MinL_kn = MinL_kn + fL*1i*(MLin_mode(:,1+(ii-1)*N:ii*N)-(MLin_mode(:,1+(ii-1)*N:ii*N))');
    end
    end
end

for ii = 1:N
    if abs(abs(Yn(ii,ii))-1) < 1e-6 
    if real(log(Yn(ii,ii))/1i) > 0
        MoutL_kp = MoutL_kp + (1-fL)*1i*(ML_mode(:,1+(ii-1)*N:ii*N)-(ML_mode(:,1+(ii-1)*N:ii*N))');
        MinR_kp = MinR_kp + fR*1i*(MRin_mode(:,1+(ii-1)*N:ii*N)-(MRin_mode(:,1+(ii-1)*N:ii*N))');
    else
        MoutL_kn = MoutL_kn + (1-fL)*1i*(ML_mode(:,1+(ii-1)*N:ii*N)-(ML_mode(:,1+(ii-1)*N:ii*N))');
        MinR_kn = MinR_kn + fR*1i*(MRin_mode(:,1+(ii-1)*N:ii*N)-(MRin_mode(:,1+(ii-1)*N:ii*N))');
    end
    end
end
        
        
Kpval = [Kpval real(Tkp)];
Knval = [Knval real(Tkn)];
Tkp
Tkn


IR_kp = trace(MoutR_kp*Gn(dim-N+1:dim,dim-N+1:dim)-MinR_kp*Gp(dim-N+1:dim,dim-N+1:dim))
IR_kn = trace(MoutR_kn*Gn(dim-N+1:dim,dim-N+1:dim)-MinR_kn*Gp(dim-N+1:dim,dim-N+1:dim))

IL_kp = trace(-MoutL_kp*Gn(1:N,1:N)+MinL_kp*Gp(1:N,1:N))
IL_kn = trace(-MoutL_kn*Gn(1:N,1:N)+MinL_kn*Gp(1:N,1:N))


end
end

% % plot((E_arr-pot_arr)./del,Kpval,'LineWidth',2);
% % hold on
% figure(1)
% plot((E_arr-pot_arr),(Kpval-Knval)./(Kpval+Knval),'LineWidth',2);
% hold on
