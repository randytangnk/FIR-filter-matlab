function [ phas, sgnl] = gmsk_modul( dat,T,M )
% GMSK modulator
% dat: Original data
% T ：One bit time
% M ：Sampling times over one bit duration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           Gaussian impulse generation          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%M = 20;     %Sampling times over one bit duration
L = 4;      %pulse duration in number bits
%T = 1;      %bit time
BT = 0.4;   %band time product
[phi, s] = pul_phas(M,L,BT,T);
Ts = T/M;
tb = [-L*T/2:Ts:L*T/2-Ts];
N = length(dat);
%%--------------------%%
%%---     figure   ---%%
%%--------------------%%
% figure(10);
% subplot(1,2,1);
% plot(tb,s);
% title('GMSK信号1比特符号波形');
% subplot(1,2,2);
% plot(tb,phi);
% title('1比特符号对应的相位');
% [m n] = size(s);
% ss = zeros(m,n);
% for(i = 1:n)
%    ss(i)=  sum(s(1:i));
% end
% ss = ss*(T/M);
% subplot(1,3,3);
% plot(ss);

phas = zeros(1,L*M);
sgnl = zeros(1,L*M);
phas_tmp = zeros(1,L*M);
sgnl_tmp = zeros(1,L*M);
% Modulation
for i = 1:N
   if dat(i) == 1
      phas_tmp = [phas(i*M+1:(i+L-1)*M) phas((i+L-1)*M)*ones(1,M)] + phi;  
     %错误逻辑：
     %phas_tmp = [phas(i*M+1:(i+L-1)*M) zeros(1,M)] + phi;
     %phas作为积分值，应该补上最后的积分值，而不应该补零！！！！！！
      sgnl_tmp = [sgnl(i*M+1:(i+L-1)*M) zeros(1,M)] + s;
   else
      phas_tmp = [phas(i*M+1:(i+L-1)*M) phas((i+L-1)*M)*ones(1,M)] - phi;
      sgnl_tmp = [sgnl(i*M+1:(i+L-1)*M) zeros(1,M)] - s;
   end
   phas = [phas(1:i*M) phas_tmp];
   sgnl = [sgnl(1:i*M) sgnl_tmp];
    
end

% %另一种求相位思路,用积分
% [m n] = size(sgnl)
% ss = zeros(m,n);
% for(i = 1:n)
%    ss(i)=  sum(sgnl(1:i));
% end
% ss = pi*ss*(T/M);
% subplot(4,1,3);
% plot(ss);
% grid on;

t = 0:Ts:(L+N)*T-Ts;
%%--------------------%%
%%---     figure   ---%%
%%--------------------%%
% figure(2);
% subplot(1,2,1)
% plot(t,sgnl);
% grid on;
% subplot(1,2,2)
% plot(t,phas);
% grid on;


end

