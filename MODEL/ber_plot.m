clc
clear all
close all


%% ================================= CRC Type
% lora_ber = load("lora_ber.mat").BER;
lora_rs_ber = load('lora_rs_ber.mat').BER;
lora_rsg_ber = load('lora_rsg_ber.mat').BER;
lora_rsg_mod_ber = load('lora_rsg_mod_ber.mat').BER;
snr = [-20:0];



figure(1)
semilogy(snr, lora_rs_ber,'-*','color','k');
hold on
semilogy(snr, lora_rsg_ber,'-s','color','r');
semilogy(snr, lora_rsg_mod_ber,'-o','color','b');
% semilogy(snr, lora_crcrs_ber,'-+','color','g');
legend('RS', 'RSG', 'RSG M')
grid
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
hold off
title('SNR');