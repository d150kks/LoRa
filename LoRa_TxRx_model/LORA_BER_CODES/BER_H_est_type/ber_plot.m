clc
clear all
close all


%% ================================= CRC Type
HestCRC_ber = load("HestCRC_ber.mat").BER;
HestNOCRC_ber = load("HestNOCRC_ber.mat").BER;
snr = load("snr.mat").snr;

figure(1)
semilogy(snr, HestCRC_ber,'-*','color','k');
hold on
semilogy(snr, HestNOCRC_ber,'-s','color','r');
legend('Hest CRC', 'Hest NOCRC')
grid on
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
hold off
title('SNR');


% lora_Hideal_F_H_ber - Ideal H, 1-freq est, 2-H est
% lora_Hest_F_H_ber - Non-Ideal H, 1-freq est, 2-H est
% lora_Hideal_H_F_ber - Ideal H, 1-H est, 2-freq est
% lora_Hest_H_F_ber - Non-Ideal H, 1-H est, 2-freq est
