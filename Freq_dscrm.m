function [ freq_out ] = Freq_dscrm( S, f0, Ts, t )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%                  ���������±�Ƶ                 %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_all = length(t);
%f1 = f0/2; 
m = 1;
%m = 1;
RI = S.*cos(m*f0*t);       % �±�Ƶͬ��·
RQ = -S.*sin(m*f0*t);      % �±�Ƶ����·
% figure(2);
% subplot(5,1,1);
% plot(t,RI);
% title('�����±�Ƶͬ��·');
% grid on;
% subplot(5,1,2);
% plot(t,RQ);
% title('�����±�Ƶ����·');
% grid on;

% Lowband pass filtering
RI = LPF((m)*f0*Ts,RI,51); %��ͨ�˲�����ֹƵ��Ϊm+0.5����Ƶ

RQ = LPF((m)*f0*Ts,RQ,51); %����·��ͬ������

% 
% subplot(5,1,3);
% plot(t,RI);
% title('��ͨ�˲���ͬ��·');
% grid on;
% subplot(5,1,4);
% plot(t,RQ);
% title('��ͨ�˲�������·');
% grid on;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%                   ���ּ�Ƶ                     %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq_out = zeros(1,N_all);
freq_out(1) = 0;
for i = 2:1:N_all
   freq_out(i) = (RQ(i)*RI(i-1)-RI(i)*RQ(i-1));%/(RI(i)^2+RQ(i)^2);
  %phase_out(i) = atan(SQ(i)/SI(i));
end
freq_out = LPF((0.5)*f0*Ts,freq_out,51);
end

