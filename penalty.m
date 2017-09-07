tetta = 0*pi;
eps_NN = 36.4;
eps_CC = 57.5;
sigma_NN = 0.3318*10^-9;
sigma_CC = 0.26*10^-9;
sigma_CN = 1/2*(sigma_CC+sigma_NN);
d = 0.1098*10^-9;
delta_C = 0.335*10^-9;
eps_CN = sqrt(eps_CC*eps_NN);
rho_C = 114*10^27;
i=0;

for r_ij=0.3:0.001:0.6
i = i+1;
dop=0;
for j = 0:1:4
    delta_r = j*delta_C;
    N1 = r_ij*10^-9+delta_r+d/2*cos(tetta);
    N2 = r_ij*10^-9+delta_r-d/2*cos(tetta);
    dop=dop+2*pi*eps_CN*rho_C*delta_C*sigma_CN^2*(2/5*(sigma_CN/N1)^10-(sigma_CN/N1)^4);
    dop=dop+2*pi*eps_CN*rho_C*delta_C*sigma_CN^2*(2/5*(sigma_CN/N2)^10-(sigma_CN/N2)^4);
end
u(i)=dop;
%u(i) = u(i)/eps_CN;
dist(i) = r_ij;
end
plot(dist,u);