clc
clear all
close all


%% ================================= Переменные
% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
rc_size = 4;
rc = (SF-rc_size);
BW = 2e6;
snr = [-6:1:10];
nIter = 100;

LORA = myLoRaClass_true(SF,BW);
LORA2 = myLoRaClass_true(SF+1,BW);

% Base = LORA.Base;
downch = LORA.downch;
chirp1 = LORA.chirp;
chirp2 = LORA2.chirp;
downch2 = LORA2.downch;

modch = [chirp1, chirp1] + chirp2;
% dp = sum(modch.*chirp2)
fourier = abs(fft(modch.*downch2));
% fourier = abs(fft(modch(1:end/2).*downch));
% fourier = abs(fft(chirp1.*downch));

figure(1)
% plot(real(chirp))
stem( 10*log10(fourier) )
    

% % rs_array = [1, 2, 3, 4];
% for rs = 1:3:4
%     LORA = myLoRaClass_true(SF,BW);
%     LORA.rs_size = rs;
%     LORA.fir_win = 1;
%     Base = LORA.Base;
%     downch = LORA.downch;
%     chirp = LORA.chirp;
%     
%     num_pre = 4;
%     
%     num_sym = 1000;
%     numinfobits = num_sym*rc; 
%     data = randi([0 1],1, numinfobits); 