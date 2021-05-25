kx_arr = -pi:0.01*pi:pi;
Ek = [];
for kx = kx_arr
    val = eig(sHk(kx));
    Ek = [Ek val];
end
sz = size(Ek);
sz = sz(1);
figure
for ii = 1:sz
    plot(kx_arr./pi,Ek(ii,:),'LineWidth',1.5)
    hold on;
end
% plot([-1 1],[E-pot E-pot]);
% hold off;