clc
clear all
close all


tic

%% ================================= Переменные
% Class initialization



% коэффициенты
SF = 8;        % коэффициент расширения спектра (от 7 до 12)
BW = 30e6;
LORA = myLoRaClass(SF, BW);

num_pre = 8;   % Число символов Преамбулы
num_check = 2; % Число символов Синхронизации
snr = [-20:2:-0];

% массивы данных
% numinfobits = 2310;          % число бит
numinfobits = 20000;          % число бит
data = randi([0, 1], 1, numinfobits); % формирование массива бит
num_sym = numinfobits/SF;



%% ================================= модуляция

[mod_chirp, check_data] = LORA.lorax_modified( data, num_sym, 0);


%% ================================= BER
tic
for n = 1:length(snr)
    n 

    [numErr, NumData] = deal(0);

    while  NumData < 2e5
        
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
save('lora_ber.mat','BER')
save('snr.mat','snr')

