clc
clear all
close all


%% ================================= CRC Type
lora_FEv1_ber = load("lora_FEv1_ber.mat").BER;
lora_FEv2_ber = load("lora_FEv2_ber.mat").BER;
lora_FEv3_ber = load("lora_FEv3_ber.mat").BER;
snr = load("snr.mat").snr;



figure(1)
semilogy(snr, lora_FEv1_ber,'-*','color','k');
hold on
semilogy(snr, lora_FEv2_ber,'-s','color','r');
semilogy(snr, lora_FEv2_ber,'-o','color','b');
legend('FEv1', 'FEv2', 'FEv3')
grid
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
hold off
title('SNR');

% FEv1 - mean peaks of fourier in CFOint
% FEv2 - sum fourier in CFOint
% FEv3 - sum fourier in CFOint and sum fourier in CFOfraq
