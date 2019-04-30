function [ freq_out ] = Freq_dscrm( S, f0, Ts, t )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%                  数字正交下变频                 %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_all = length(t);
%f1 = f0/2; 
m = 1;
%m = 1;
RI = S.*cos(m*f0*t);       % 下变频同相路
RQ = -S.*sin(m*f0*t);      % 下变频正交路
% figure(2);
% subplot(5,1,1);
% plot(t,RI);
% title('数字下变频同相路');
% grid on;
% subplot(5,1,2);
% plot(t,RQ);
% title('数字下变频正交路');
% grid on;

% Lowband pass filtering
RI = LPF((m)*f0*Ts,RI,51); %低通滤波，截止频率为m+0.5倍载频

RQ = LPF((m)*f0*Ts,RQ,51); %正交路作同样处理

% 
% subplot(5,1,3);
% plot(t,RI);
% title('低通滤波后同相路');
% grid on;
% subplot(5,1,4);
% plot(t,RQ);
% title('低通滤波后正交路');
% grid on;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%                   数字鉴频                     %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq_out = zeros(1,N_all);
freq_out(1) = 0;
for i = 2:1:N_all
   freq_out(i) = (RQ(i)*RI(i-1)-RI(i)*RQ(i-1));%/(RI(i)^2+RQ(i)^2);
  %phase_out(i) = atan(SQ(i)/SI(i));
end
freq_out = LPF((0.5)*f0*Ts,freq_out,51);
end

