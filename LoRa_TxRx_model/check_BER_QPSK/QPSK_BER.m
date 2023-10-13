clc
clear all
close all


tic

%% ================================= Переменные

% коэффициенты
M = 4;
snr = [-20:2:20];

numinfobits = 1e4;          % число бит
data = randi([0, 1], 1, numinfobits); % формирование массива бит


%% ================================= модуляция
mod_sig = qammod(data.', M, 'InputType','bit').';

%% ================================= BER
tic
for n = 1:length(snr)
    n 

    [numErr, NumData] = deal(0);

    while  NumData < 1e5
        
        % АБГШ 
        rxSig = awgn(mod_sig,snr(n),'measured');

        % демодуляция
        demod_sig = qamdemod(rxSig.', M, 'OutputType','bit').';

        % подсчет БЕР с учетом задержки
        err = sum(demod_sig~=data);

        % Increment the error and bit counters
        numErr = numErr + err;        
        NumData = NumData + numinfobits;
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
save('qpsk_ber.mat','BER')
save('snr.mat','snr')

