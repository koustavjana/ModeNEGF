M11 = -eye(N);
M12 = inv(H21)*(E*eye(N)-H22);
M21 = -1*inv(H12)*(E*eye(N)-H11);
M22 = inv(H12)*(E*eye(N)-H11)*inv(H21)*(E*eye(N)-H22) - eye(N);

M = [M11 M12; M21 M22];
[V,D] = eig(M);
D = diag(D);

kD = log(D)/1i;

Up = [];
Un = [];
for ii=1:2*N
    if abs(imag(kD(ii))) > 1e-10
        if imag(kD(ii))>0
            Up=[Up V(:,ii)];
        end
        if imag(kD(ii))<0
            Un=[Un V(:,ii)];
        end
    end
    if abs(imag(kD(ii))) <= 1e-10
        vel = V(:,ii)'*vk(kD(ii))*V(:,ii);
        if vel > 0
            Up=[Up V(:,ii)];
        end
        if vel <= 0
            Un=[Un V(:,ii)];
        end
    end
end
Up1 = Up(1:N,:);
Up2 = Up(N+1:2*N,:);
Un1 = Un(1:N,:);
Un2 = Un(N+1:2*N,:);
Yp = diag(Up'*M*Up);
Yn = diag(Un'*M*Un);
Yp = diag(Yp);
Yn = diag(Yn);
Fp1 = Up1*Yp*inv(Up1);
Fp2 = Up2*Yp*inv(Up2);
Fn1 = Un1*Yn*inv(Un1);
Fn2 = Un2*Yn*inv(Un2);

Ukp = [];
Ukn = [];
for ii = 1:2*N
   if real(kD(ii))~=pi && real(kD(ii))~=0
    if real(kD(ii))>0
        Ukp=[Ukp V(:,ii)];
    end
    if real(kD(ii))<=0
        Ukn=[Ukn V(:,ii)];
    end
   end
   if real(kD(ii))==pi || real(kD(ii))==0
    if imag(kD(ii))>0
        Ukp=[Ukp V(:,ii)];
    end
    if imag(kD(ii))<=0
        Ukn=[Ukn V(:,ii)];
    end
   end
end
Ykp = diag(Ukp'*M*Ukp);
Ykn = diag(Ukn'*M*Ukn);
Ykp = diag(Ykp);
Ykn = diag(Ykn);
Ukp1 = Ukp(1:N,:);
Ukp2 = Ukp(N+1:2*N,:);
Ukn1 = Ukn(1:N,:);
Ukn2 = Ukn(N+1:2*N,:);
