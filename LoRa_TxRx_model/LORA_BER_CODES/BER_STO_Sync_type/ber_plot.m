clc
clear all
close all


%% ================================= CRC Type
lora_STOv1_ber = load("lora_STOv1_ber.mat").BER;
lora_STOv2_ber = load("lora_STOv2_ber.mat").BER;
lora_STOv1_H_ber = load("lora_STOv1_H_ber.mat").BER;
lora_STOv2_H_ber = load("lora_STOv2_H_ber.mat").BER;
snr = load("snr.mat").snr;

figure(1)
semilogy(snr, lora_STOv1_ber,'-*','color','k');
hold on
semilogy(snr, lora_STOv2_ber,'-s','color','r');
semilogy(snr, lora_STOv1_H_ber,'-*','color','b');
semilogy(snr, lora_STOv2_H_ber,'-s','color','g');
legend('STOv1', 'STOv2', 'STOv1 H', 'STOv2 H')
grid on
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
hold off
title('SNR');



% subplot(121)
% semilogy(snr, lora_STOv1_ber,'-*','color','k');
% hold on
% semilogy(snr, lora_STOv2_ber,'-s','color','r');
% legend('STOv1', 'STOv2')
% grid
% xlabel('SNR (dB)')
% ylabel('Bit Error Rate')
% hold off
% title('SNR');
% 
% subplot(122)
% semilogy(snr, lora_STOv1_H_ber,'-*','color','k');
% hold on
% semilogy(snr, lora_STOv2_H_ber,'-s','color','r');
% legend('STOv1 H', 'STOv2 H')
% grid
% xlabel('SNR (dB)')
% ylabel('Bit Error Rate')
% hold off
% title('SNR');

% STOv1 - 1 downch x 1 preamble
% STOv2 - 1 downch x 8 preamble
% STOv1_H - 1 downch x 1 preamble with channel H
% STOv2_H - 1 downch x 8 preamble with channel H
