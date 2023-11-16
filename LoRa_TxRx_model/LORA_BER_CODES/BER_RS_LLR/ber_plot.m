clc
clear all
close all


%% ================================= CRC Type
lora_rsg_ber = load("lora_rsg_ber.mat").BER;
lora_ber = load("lora_ber.mat").BER;
snr = [-16:0];



figure(1)
semilogy(snr, lora_rsg_ber,'-*','color','k');
hold on
semilogy(snr, lora_ber,'-s','color','r');
legend('RSG', 'Classic')
grid
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
hold off
title('SNR');
