clc
clear all
close all


%% ================================= CRC in Noise channel
lora_ber = load("lora_crc_ber.mat").BER;
lora_crc_ldpc_ber = load("lora_crc_ldpc_ber_1944.mat").BER;
snr = load("snr_crc.mat").snr;



figure(1)
semilogy(snr, lora_ber,'-*','color','k');
hold on
semilogy(snr, lora_crc_ldpc_ber,'-s','color','r');
legend('Without LDPC', 'With LDPC')
grid
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
hold off
title('SNR');
