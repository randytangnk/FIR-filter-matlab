function [ phi,s ] = pul_phas( T_over_Ts, L, BT, T )
% GMSK pulse generator function
% T_over_Ts : Sample times over one bit duration
% L : Pulse duration number of bit time 
% BT : Band time pruduct
% T: One bit time

Ts = T/T_over_Ts;
Nphi = ceil(T_over_Ts*L);
Nts = ceil(T_over_Ts);
phi = zeros(1, Nphi);
s = zeros(1,Nphi);

sigma = sqrt(log(2))/2/pi/(BT/T);

t = [-L*T/2:Ts:L*T/2];
t = t(1:Nphi);
ta = t + T/2;
tb = t - T/2;

Qta = 0.5*(ones(1,Nphi) + erf(ta/sigma/(sqrt(2))));
Qtb = 0.5*(ones(1,Nphi) + erf(tb/sigma/(sqrt(2))));
expta = exp(-0.5*((ta/sigma).^2))/sqrt(2*pi)*sigma;
exptb = exp(-0.5*((tb/sigma).^2))/sqrt(2*pi)*sigma;

phi = pi/T/2*(ta.*Qta+expta - tb.*Qtb - exptb);

s = 1/2/T*(Qta - Qtb);

% figure(20);
% subplot(1,2,1)
% plot(t,s);
% grid on;
% subplot(1,2,2)
% plot(t,phi);
% grid on;

end

