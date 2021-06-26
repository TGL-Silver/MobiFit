function [P1]=Ceullar(F,data)
Fs = F;
T = 1/Fs;
L = length(data);
h = hamming(L, 'periodic');
data=data.*h;
Y = fft(data);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% global s;
% 
% s=officehuo;


% A=((sum(F(2,CellularS(count):(CellularPd-1))).^2)/((CellularPd-1)-CellularS(count)))
% B=((sum(F(2,CellularE(count-1):CellularS(count))).^2)/(CellularS(count)-CellularE(count-1)))
