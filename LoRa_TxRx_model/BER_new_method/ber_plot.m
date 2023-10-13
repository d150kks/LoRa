clc
clear all
close all

%% ================================= Переменные
% lora_ber = load("lora_ber.mat").BER;
% lora_crc_ber = load("lora_crc_ber.mat").BER;
% qpsk_ber = load("qpsk_ber.mat").BER;
% snr = load("lora_snr.mat").snr;


% lora_ber = load("lora_ber.mat").BER;
lora_crc_pluto_ber = load("lora_crc_pluto_ber.mat").BER;
% qpsk_ber = load("qpsk_ber.mat").BER;
snr = load("lora_crc_pluto_snr.mat").tx_power;


figure(1)
% semilogy(snr,lora_ber,'-*','color','k');
% hold on
semilogy(snr,lora_crc_pluto_ber,'-s','color','r');
% semilogy(snr,qpsk_ber,'-o','color','b');
% legend('lora', 'qpsk')
grid
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
hold off
title('SNR');