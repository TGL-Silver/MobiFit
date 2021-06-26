function [S,W,T,tax]=NewFreCutSTFT(sNorm,nsc,nov,nff,fs)


h = hamming(nsc, 'periodic');   % Hamming weight function
L = length(sNorm);              % Length of Signal
nst = nsc-nov;                  % Number of STep per colum
% fix是指靠近0取整
ncl = fix( (L-nsc)/nst ) + 1;   % Number of CoLum
%ncl = floor( (L-nsc)/nst ) + 1;
% nextpow2是指靠的最近的2的指数，用在需要数组长度为2的指数的情况
%nff = 2^nextpow2(nsc);         % Number of Fast Fourier Transformation
% nff = max(256,2^nextpow2(nsc));
Ps = zeros(nsc,ncl);
for n = 1:ncl                   % Ps means Processed Signal
    %         size(signal( (n-1)*nst+1 : (n-1)*nst+nsc ))
    %         size(h)
    
    %%  设置与不设置海明窗区别
    Ps(:,n) = sNorm( (n-1)*nst+1 : (n-1)*nst+nsc ).*h';
    %
%     Ps(:,n) = sNorm( (n-1)*nst+1 : (n-1)*nst+nsc );
    
end                             % Ps is a matrix

% Short-time Fourier transform
STFT0 = fft(Ps,nff);
% Turn into DFT in dB
STFT1 = abs(STFT0/nff);
STFT2 = STFT1(1:nff/2+1,:);             % Symmetric(对称的) results of FFT
STFT2(2:end-1,:) = 2*STFT2(2:end-1,:);  % Add the value of the other half
% STFT3 = 20*log10(STFT2);                % Turn sound pressure into dB level

% Axis Generating
fax = fs*(0:(nff/2))/nff;                           % Frequency axis setting
tax = ( .5*nsc : nst : nst*(ncl-1)+.5*nsc ) / fs;   % Time axis generating
% meshgrid 生成网格矩阵
[ffax,ttax] = meshgrid(tax,fax);                    % Generating grid of figure
% Output
W = ffax;
T = ttax;
% Snormalize=normalize(STFT2(2:end,:),0,1);
% S = [STFT2(1,:);Snormalize];
S=STFT2;
figure;
my_pcolor(W(2:20,:),T(2:20,:),S(2:20,:));
title('STFT','FontSize',18,'fontname','Times New Roman');
xlabel('Time(seconds)','FontSize',18,'Fontname', 'Times New Roman');
ylabel(' Frequency ','FontSize',18,'Fontname', 'Times New Roman');
set(gca,'FontSize',18)



end
