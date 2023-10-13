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


%% ================================= BER
tic

for n = 1:length(snr)
    n 

    [numErr, NumData] = deal(0);

    for iter = 1:nIter
        
        % АБГШ 
        rxSig = awgn(mod_chirp,snr(n),'measured');

        % демодуляция
        [sv_decode, sv, fourier] = LORA.delorax_modified(rxSig, num_sym);
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
save('lora_ber.mat','BER')
save('snr.mat','snr')

