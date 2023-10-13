clc
clear all
close all


tic

%% ================================= Переменные
% Class initialization



% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
BW = 30e6;

LORA = myLoRaClass(SF, BW);
Base = LORA.Base;
downch = LORA.downch;

snr = [-16:1:0];
nIter = 10;

% массивы данных
nbits = 11200; 
num_sym = nbits/SF;
data = randi([0 1],1, nbits); 

%% ================================= Rate matching
[dataRM, numcodebitsRM, num_sym, zeros2end, flagRM] = LORA.RM(data);


%% ================================= CRC coding
[data_code] = LORA.codeCRC(dataRM, num_sym);

%% ================================= модуляция
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified( data_code, num_sym, 1);
num_pre = 8;
tx_preamble = repmat(LORA.chirp,1,num_pre);
tx_chirp = [downch, tx_preamble, mod_chirp]; % DOWN, UP, PAYLOAD
% tx_chirp = [tx_preamble, mod_chirp]; % DOWN, UP, PAYLOAD
tx_length = length(tx_chirp);
sig_length = tx_length-Base;

%% ================================= BER
tic

h11 = zeros(1, Base);
h11(1) = 1;
h11(2) = 0.5;
h11(5) = 0.3;
H11 = fft(h11);

% tx_chirp = ifft( fft(tx_chirp).*H11 );
% tx_chirp = tx_chirp;
for i=1:tx_length/Base
    tx_sig_fft = fft(tx_chirp(i*Base-Base+1:Base*i));

    tx_chirp(i*Base-Base+1:Base*i) = ifft(tx_sig_fft.*H11);
end

fps = BW/Base;
freq_shift = -fps*1.5;
dphi=freq_shift*2*pi*(1/BW);% сдвиг

% вводим частотный сдвиг
for j=1:length(tx_chirp)
    tx_chirp(j)=tx_chirp(j)*exp(1i*dphi*j);
end

% BER Calculation
for n = 1:length(snr)
    n 

    [numErr, NumData] = deal(0);

    for iter = 1:nIter
        
        % АБГШ 
        delay = randi([1e2,1e3]); 
        rxSig_noncor = awgn([zeros(1,delay), tx_chirp,zeros(1,500)],snr(n),'measured');

        % Correlation
        [cor,lags] = xcorr(rxSig_noncor, tx_preamble);
        cor(round(end/2)+10:end) = cor(round(end/2)+10:end).*gausswin(length(cor(round(end/2)+10:end))).';
        [~, max_idx] = max(abs(cor));
        start = lags(max_idx);
        if(start>delay)
            start = delay-2;
        end
        rxSig = rxSig_noncor(abs(start)+1:end);
%         rxSig = awgn([tx_chirp,zeros(1,500)],snr(n),'measured');


        % Frequency correction
        [freq_data, corrected_signal, rx_preamble] = LORA.LORA_FREQ_ESTIM_v2(rxSig, num_pre, sig_length);

        % демодуляция
        [soft_bits, hard_bits, sv, sv_cor, fourier] = LORA.DELORAX_CRC(corrected_signal, num_sym, tx_preamble, rx_preamble);
        
        % CRC decoding
        [data_decode] = LORA.decodeCRC(hard_bits, num_sym, zeros2end, flagRM);

        % подсчет БЕР с учетом задержки
        err = sum(data_decode~=data);

        % Increment the error and bit counters
        numErr = numErr + err;
        NumData = NumData + nbits;
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

% return
% 
% save('lora_STOv1_ber.mat','BER')
% save('lora_STOv2_ber.mat','BER')
% save('lora_STOv1_H_ber.mat','BER')
% save('lora_STOv2_H_ber.mat','BER')
save('snr.mat','snr')

