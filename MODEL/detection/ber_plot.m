clc
clear all
close all

snr = [-16:0];
ber_fft = load('BER_FFTSYNCH.mat').BER;
ber_snr = load('BER_snrSYNCH.mat').BER;


figure(1)
semilogy(snr,ber_fft,'-*','color','k');
hold on
semilogy(snr,ber_snr,'-*','color','b');
grid
legend('FFT прием', 'Корреляционный прием')
xlabel('ОСШ (дБ)')
ylabel('Вероятность битовой ошибки')
hold off