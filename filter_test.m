clc;
clear;

M = 25;     %Sample 20 times during one bit time
L = 4;      %Gaussian pulse duration for one bit 
T = 1;      %One bit time supposed to be 1
Ts = T/M;   %Sampling period
% Generate random binary data
N = 20; % 20 bit long data

%dat0 = [0 0 1 1 1 0 1 1 0 0 1 1 0 0 0 0 1 1 1 1];
dat0 = repmat([0 0 1 1],1,N/4);  %dat0重复的01序列
%dat0 = rand(1,N);
%dat0 = repmat([0 1],1,N/2);
%dat1 = rand(1,N);
dat1 = repmat([0 0],1,N/2);
dat2 = rand(1,N);
%dat2 = repmat([0 1],1,N/2);
dat3 = rand(1,N);
%dat3 = repmat([0 1],1,N/2);

% BPSK modulation
for i = 1:N
   if dat0(i) > 0.5
       dat0(i) = 1;
   else
       dat0(i) = -1;
   end
   if dat1(i) > 0.5
       dat1(i) = 1;
   else
       dat1(i) = -1;
   end
   if dat2(i) > 0.5
       dat2(i) = 1;
   else
       dat2(i) = -1;
   end
   if dat3(i) > 0.5
       dat3(i) = 1;
   else
       dat3(i) = -1;
   end

end

t = 0:Ts:(L+N)*T-Ts;   %time axis ，totally (L+N)*T/Ts = （4 + 20）*1/(1/20) = 480 samples
f0 = 2*pi*2/T;       % Center frequency,1/T = 9.6kHz 

[phas0 sgnl0] = gmsk_modul(dat0,T,M); %GMSK调制产生相位和信号波形
df0 = 0;
f0 = 2*pi*25/9.6/T;
S0 = cos((1+df0)*f0*t+phas0+0.12*pi); 
N_all = length(phas0); %480 total samples

a1 = 1;  %假设接收到的其他路信号衰减0.5，时延60个采样点
d1 = 0;
df1 = 0;
f1 = 2*pi*50/9.6/T;
[phas1 sgnl1] = gmsk_modul(dat1,T,M);
S1 = cos((1+df1)*f1*t+phas0);
S1 = a1*[zeros(1,d1) S1(1:(N_all-d1))];

a2 = 1;  %衰减0.3，时延79个采样点
d2 = 0;
df2 = 0;
f2 = 2*pi*75/9.6/T;
[phas2 sgnl2] = gmsk_modul(dat2,T,M);
S2 = cos((1+df2)*f2*t+phas0);
S2 = a2*[zeros(1,d2) S2(1:(N_all-d2))];
S = S0 +S1 +S2;
%test filter
SL = LPF((f0+f1)/2*Ts,S,51);
SB = BPF((f0+f1)/2*Ts,(f1+f2)/2*Ts,S,51);
SH = HPF((f1+f2)/2*Ts,S,51);
SS = BSF((f0+f1)/2*Ts,(f1+f2)/2*Ts,S,51);
%FD
SLF = Freq_dscrm( SL, f0, Ts, t );
SBF = Freq_dscrm( SB, f1, Ts, t );
SHF = Freq_dscrm( SH, f2, Ts, t );

figure(1)
subplot(2,1,1);
plot(S);
title('三路信号时域');
subplot(2,1,2);
plot(abs(fft(S)));
title('三路信号频域');

figure(2)
subplot(4,2,1);
plot(SL);
title('低通滤波后时域');
subplot(4,2,2);
plot(abs(fft(SL)));
title('低通滤波后频域');
subplot(4,2,3);
plot(SB);
title('带通滤波后时域');
subplot(4,2,4);
plot(abs(fft(SB)));
title('带通滤波后频域');
subplot(4,2,5);
plot(SH);
title('高通滤波后时域');
subplot(4,2,6);
plot(abs(fft(SH)));
title('高通滤波后频域');
subplot(4,2,7);
plot(SS);
title('带阻滤波后时域');
subplot(4,2,8);
plot(abs(fft(SS)));
title('带阻滤波后频域');

figure(3);
subplot(3,1,1);
plot(SLF);
title('低通滤波后鉴频');
subplot(3,1,2);
plot(SBF);
title('带通滤波后鉴频');
subplot(3,1,3);
plot(SHF);
title('高通滤波后鉴频');






