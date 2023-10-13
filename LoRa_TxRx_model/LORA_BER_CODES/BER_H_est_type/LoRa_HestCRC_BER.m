clc
clear all
close all


tic

%% ================================= Переменные
% Class initialization
% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
BW = 125e3;

LORA = myLoRaClass(SF, BW);
Base = LORA.Base;
downch = LORA.downch;
chirp = LORA.chirp;
Ts = LORA.Ts;

snr = [-16:1:0];
nIter = 100;

% массивы данных
nbits = 1200; 
num_sym = nbits/(SF-4);
data = randi([0 1],1, nbits); 

%% ================================= Rate matching
[dataRM, numcodebitsRM, num_sym, zeros2end, flagRM] = LORA.RM(data);

%% ================================= CRC coding
[data_code] = LORA.codeCRC(dataRM, num_sym);

%% ================================= модуляция
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified( data_code, num_sym, 1);
num_pre = 8;
tx_preamble = repmat(LORA.chirp,1,num_pre);
tx_downch = repmat(downch,1,num_pre);
tx_chirp = [tx_downch, tx_preamble, mod_chirp]; % DOWN, UP, PAYLOAD
% tx_chirp = [downch, tx_preamble, mod_chirp]; % DOWN, UP, PAYLOAD
tx_length = length(tx_chirp);

%% ================================= BER
tic

h11 = zeros(1, Base);
h11(1) = 1;
h11(2) = 0.9;
h11(3) = 0.8;
h11(5) = 0.8;
h11(7) = 0.7;
h11(9) = 0.5;
h11(11) = 0.45;
h11(13) = 0.3;
h11(15) = 0.2;
h11(16) = 0.1;
H11 = fft(h11);
% H11 = ones(1, Base);

for i=1:tx_length/Base
    tx_sig_fft = fft(tx_chirp(i*Base-Base+1:Base*i));
    tx_chirp_h(i*Base-Base+1:Base*i) = ifft(tx_sig_fft.*H11);
end
% tx_chirp = tx_chirp;

fps = BW/Base;
freq_shift = fps*1.5;
dphi=freq_shift*2*pi*(1/BW);% сдвиг

% вводим частотный сдвиг
for j=1:length(tx_chirp_h)
    tx_chirp_hf(j)=tx_chirp_h(j)*exp(1i*dphi*j);
end

% BER Calculation
% snr=20;
for n = 1:length(snr)

    [numErr, NumData] = deal(0);

    for iter = 1:nIter
        fprintf( 'iter num: %.d\n', n); 
        
        % АБГШ 
        delay = 2;
        rxSig = awgn([zeros(1,delay), tx_chirp_hf,zeros(1,100)],snr(n),'measured');
        rxSig = rxSig(1:tx_length);

        [freq_data, corrected_signal, rx_preamble] = LORA.LORA_FREQ_ESTIM_v3(rxSig, num_pre);
        fprintf( 'Err:     %.2f\n', freq_shift-(freq_data{2}+freq_data{3}+freq_data{4}))
        fprintf( 'STO:     %.2f\n', freq_data{1});


H11hat = 0;
for i=1:num_pre
    H11hat = H11hat + fft(rx_preamble(i*Base-Base+1:i*Base))./fft(chirp);
end
H11hat = H11hat/(num_pre);

for i=1:length(corrected_signal)/Base
    symrxfft = fft(corrected_signal(i*Base-Base+1:Base*i));
    corrected_signal(i*Base-Base+1:Base*i) = ifft( symrxfft./H11hat );
%     corrected_signal(i*Base-Base+1:Base*i) = ifft( symrxfft./H11 );
end
% nvar = std( abs(mod_chirp-corrected_signal).^2);
% fprintf( 'nvar:     %.2f\n', nvar);

% figure(1)
% subplot(211); hold on
% plot(real(H11))
% plot(real(H11hat))
% 
% subplot(212); hold on
% plot(imag(H11))
% plot(imag(H11hat))
% 
% figure(2); hold on
% plot(abs(H11))
% plot(abs(H11hat))
% % 
% % figure(3); hold on
% % plot(real(corrected_signal))
% % plot(real(mod_chirp))
% return

        % демодуляция
        [soft_bits, hard_bits, sv, sv_cor, fourier] = LORA.DELORAX_CRC(corrected_signal, num_sym, tx_preamble, rx_preamble);

% figure(4)
% (stem(fourier))
% return

        % CRC decoding
        [data_decode] = LORA.decodeCRC(hard_bits, num_sym, zeros2end, flagRM);

        % подсчет БЕР с учетом задержки
        err = sum(data_decode~=data);

        % Increment the error and bit counters
        numErr = numErr + err;
        NumData = NumData + nbits;
        clc;
    end


    % Estimate the BER for both methods
    BER(n) = numErr/NumData;
end
toc

%%
figure(1)
semilogy(snr,BER,'-*','color','k');
grid
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
hold off
title('SNR');

% save('HestCRC_ber.mat','BER')
% save('snr.mat','snr')

