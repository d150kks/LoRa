clc
clear all
close all


%% ================================= CRC Type
lora_ber = load("lora_ber.mat").BER;
lora_crc_ber = load("lora_crc_ber.mat").BER;
% lora_rs_ber = load('lora_rs_ber.mat').BER;
lora_rs_ber = load('lora_rs_ber2.mat').BER;
lora_crcrs_ber = load('lora_crcrs_ber.mat').BER;
snr = load("snr_crc.mat").snr;



figure(1)
semilogy(snr, lora_ber,'-*','color','k');
hold on
semilogy(snr, lora_crc_ber,'-s','color','r');
semilogy(snr, lora_rs_ber,'-o','color','b');
semilogy(snr, lora_crcrs_ber,'-+','color','g');
legend('Without CRC', 'With CRC', 'With RS', 'With CRCRS')
grid
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
hold off
title('SNR');
