clc
clear all
close all


% upch*upch/upch*downch
% eight consequtive symbols
% large window

%% ================================= Переменные
% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
rc_size = 4;
rc = (SF-rc_size);
BW = 125e3;
nIter = 10000;
snr = 0;

LORA = myLoRaClass_true(SF,BW);
Base = LORA.Base;
downch = LORA.downch;
chirp = LORA.chirp;


%% ================================= ssk
M = 2^SF;
n = 0:M-1;
m = 13;
x0 = exp(1i*pi*(n).^2/M);
xm = exp(1i*pi*(n+m).^2/M);

% sm = exp(1i*pi*(n+4).^2/M);
x0m = x0.*x0(14).*exp(1i*2*pi*13*n/M);


fft(1:4)
(1:4).*exp(-1i*2*pi*(1:4)/4)

% Channel
% w = 0.7071*(normrnd(0,1,[1,M])+1i*normrnd(0,1,[1,M]));
w = 0;
h = 0.7071*(normrnd(0,1,[1,M])+1i*normrnd(0,1,[1,M]));
h=1;
% sm = exp(1i*pi*(n+randi(M)).^2/M);
sm = exp(1i*pi*(n+14).^2/M);
r = h.*sm + w;

R1 = fft(r.*conj(x0));
R2 = r.*x0.*exp(1i*2*pi*0*n/M);
% for i=1:M
%     mest(i) = max( sum(abs(r.*conj(circshift(x0,i)) )));
% %     mest(i) = max( real(conj(x0(i)).*R1(i)*conj(h)) );
% %     mest2(i) = max( real(x0(i).*R2(i)*conj(h)) );
% end

% A1 = abs(sum( (r.*conj(sm)).^2 ))
% A2 = abs(sum( (r.*conj(x0)).^2 ))

figure(1); hold on
plot( abs(R1))
plot( abs(R2))
% % plot( abs(R1))
% plot( abs(A1))
% plot( real(xm))
% plot( real(r))







