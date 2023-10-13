clc
clear all
close all


%% ================================= CRC Type
lora_crc_ber = load("lora_crc_ber.mat").BER;
lora_crc_type2_ber = load("lora_crc_type2_ber.mat").BER;
snr = load("snr_crc.mat").snr;



figure(1)
semilogy(snr, lora_crc_ber,'-*','color','k');
hold on
semilogy(snr,lora_crc_type2_ber,'-s','color','r');
legend('Sort by amp', 'Sort by amp and position')
grid
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
hold off
title('SNR');
