clc
clear all
close all


%% ================================= CRC Type
lora_ber = load("lora_ber.mat").BER;
lora_crc_ber = load("lora_crc_ber.mat").BER;
snr = load("snr_crc.mat").snr;



figure(1)
semilogy(snr, lora_ber,'-*','color','k');
hold on
semilogy(snr, lora_crc_ber,'-s','color','r');
legend('Without CRC', 'With CRC')
grid
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
hold off
title('SNR');
