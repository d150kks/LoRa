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
nIter = 100;

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

h11 = zeros(1, length(tx_chirp));
h11(1) = 1;
h11(2) = 0.5;
h11(5) = 0.3;
H11 = fft(h11);

% tx_chirp = ifft( fft(tx_chirp).*H11 );
% tx_chirp = tx_chirp;

fps = BW/Base;
freq_shift = fps*1.5;
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
        delay = randi(1);
        rxSig = awgn([tx_chirp,zeros(1,100)],snr(n),'measured');
        rx_preamble = rxSig(Base+1:Base*(num_pre+1));
        corrected_signal = rxSig(Base*(num_pre+1)+1:end);

        % Frequency correction
        [freq_data, corrected_signal, rx_preamble] = LORA.LORA_FREQ_ESTIM(rxSig, num_pre, sig_length);

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

%     STOint = freq_data{1};
%     est1 = freq_data{2};
%     est2 = freq_data{3};
%     est3 = freq_data{4};
%     est4 = freq_data{4};
%     freq_shift_est = est1+est2+est3+est4;
% 
%     fprintf('Num Err: %d\n', numErr)
%     fprintf('STOint:  %d\n\n', STOint)
% 
%     fprintf('Offset:    %.2f\n', freq_shift)
%     fprintf('Est 1:    %.2f\n', est1)
%     fprintf('Est 2:    %.2f\n', est2)
%     fprintf('Est 3:    %.2f\n', est3)
%     fprintf('Est 4:    %.2f\n', est4)
%     fprintf('Offset Est:  %.2f\n', freq_shift_est)
%     fprintf('Offset Err:  %.2f\n\n', abs(freq_shift_est-freq_shift))

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
save('lora_FEv1_ber.mat','BER')
save('snr.mat','snr')

