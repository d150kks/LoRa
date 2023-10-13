clc
clear all
close all


%% ================================= CRC Type
STOCFOFine_old_ber = load("STOCFOFine_old_ber.mat").BER;
STOCFOFine_ber = load("STOCFOFine_ber.mat").BER;
STOCFOFine_old_Hest_ber = load("STOCFOFine_old_Hest_ber.mat").BER;
STOCFOFine_Hest_ber = load("STOCFOFine_Hest_ber.mat").BER;
snr = load("snr.mat").snr;

STOCFOFine_Hest_NLMS_ber = load("STOCFOFine_Hest_NLMS_ber.mat").BER;


figure(1)
semilogy(snr, STOCFOFine_old_ber,'-*','color','k');
hold on
semilogy(snr, STOCFOFine_ber,'-s','color','r');
semilogy(snr, STOCFOFine_old_Hest_ber,'-o','color','k');
semilogy(snr, STOCFOFine_Hest_ber,'-+','color','r');
semilogy(snr, STOCFOFine_Hest_NLMS_ber,'-+','color','b');
legend('Fest old', 'Fest new', 'Fest old + Hest', 'Fest new + Hest', 'Fest new + Hest NLMS', 'Location','southwest')
grid
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
hold off
title('SNR');

% figure(2)
% semilogy(snr, STOCFOFine_old_Hest_ber,'-*','color','k');
% hold on
% semilogy(snr, STOCFOFine_Hest_ber,'-s','color','r');
% legend('Fest old without H', 'Fest new without H')
% grid
% xlabel('SNR (dB)')
% ylabel('Bit Error Rate')
% hold off
% title('SNR');

% STOCFOFine_old_ber - CFOfraq estimated in the end of the algorithm
% STOCFOFine_ber - CFOfraq estimated at the start of the algorithm
