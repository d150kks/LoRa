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



%% ================================= модуляция
[mod_chirp, check_data] = LORA.lorax_modified( data, num_sym, 0);

num_pre = 8;
tx_preamble = repmat(LORA.chirp,1,num_pre);
tx_chirp = [tx_preamble, mod_chirp]; % DOWN, UP, PAYLOAD

%% ================================= BER
tic

h11 = zeros(1, length(tx_chirp));
h11(1) = 1;
h11(2) = 0.5;
h11(5) = 0.3;
H11 = fft(h11);

tx_chirp = ifft( fft(tx_chirp).*H11 );

fps = BW/Base;
freq_shift = fps*0.5;
dphi=freq_shift*2*pi*(1/BW);% сдвиг

% вводим частотный сдвиг
for j=1:length(tx_chirp)
    tx_chirp(j)=tx_chirp(j)*exp(1i*dphi*j);
end

for n = 1:length(snr)
    n 

    [numErr, NumData] = deal(0);

    for iter = 1:nIter
        
        % АБГШ 
        rxSig = awgn(tx_chirp,snr(n),'measured');
        rx_preamble = rxSig(1:Base*num_pre);
        rx_sig = rxSig(Base*num_pre+1:end);

        % демодуляция
        [rx_sig_correct] = LORA_RLS(14, tx_preamble, rx_preamble, rx_sig);
        [sv_decode, sv, fourier] = LORA.delorax_modified(rx_sig_correct, num_sym);
        dbits = de2bi(sv,SF).';
        dbits = dbits(:).';

        % подсчет БЕР с учетом задержки
        err = sum(dbits~=data);

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
save('lora_rls_ber.mat','BER')
save('snr.mat','snr')

