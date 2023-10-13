clc
clear all
close all


%% ================================= CRC Search Type
lora_crc_type2_ber = load("lora_crc_type2_ber.mat").BER;
lora_crc_type2_ber_search2 = load("lora_crc_type2_ber_search2.mat").BER;
snr = load("snr_crc.mat").snr;



figure(1)
semilogy(snr, lora_crc_type2_ber,'-*','color','k');
hold on
semilogy(snr,lora_crc_type2_ber_search2,'-s','color','r');
legend('Search Local', 'Search Global')
grid
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
hold off
title('SNR');
