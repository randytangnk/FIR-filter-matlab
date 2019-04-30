function [ Y ] = BPF( omega_l, omega_h, X, N )
% BPFband pass filter using Blackman Window
% 'omega_l' is the lower cutoff radian frequency
% 'omega_h' is the higher cutoff radian frequency
% 'X' is the input sequence to be filtered.
% N is the length of FIR filter
% 'Y' is the filtered output.
% Default using Hamming window

M = length(X);
L = M+N-1;
n = 1:1:N;

h = (sin(omega_h*(n-(N+1)/2)) - sin(omega_l*(n-(N+1)/2)))/pi./(n-(N+1)/2); % Move the symetric impulse to the positive part of the time axis 
h((N+1)/2) = (omega_h-omega_l)/pi; %让其不为NAN

% Window function 
alpha = 0.5;
w = alpha + (1-alpha)*cos(2*pi/N*(n-(N+1)/2)); %Hamming Window
h = h.*w;
% Normalize the cofficient
mod = 0;
for i = 1:length(h)
   mod = mod + h(i).*h(i); 
end
mod = sqrt(mod);
h = h/mod;

% figure;
% subplot(2,1,1);
% plot(h)
% subplot(2,1,2);
% plot(abs(fft([h zeros(1,205)])));

% 将X逆序存储，后面的值表示较早时间产生，方便后面处理
for i = 1:ceil(M/2)
    temp = X(i);
    X(i) = X(M-i+1);
    X(M-i+1) = temp;
end
%Y = conv(X,h);
%Calculte convolution
Y = zeros(1,L);
X = [zeros(1,N-1) X zeros(1,N-1)];
LX = length(X);

for i = 1:L
    Y(i) = X(LX-i+2-N:LX-i+1)*h';
end
Y = Y(N-1/2:N-1/2+M-1); % 去掉卷积结果的头尾，保持和原序列等长
end

