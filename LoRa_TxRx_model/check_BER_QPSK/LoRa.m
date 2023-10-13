clc
clear all
close all


%% ================================= Переменные

% коэффициенты
SF = 8;        % коэффициент расширения спектра (от 7 до 12)
Base = 2^SF;   % База сигнала
num_sym = 100;% Число передаваемых символов
num_pre = 8;   % Число символов Преамбулы
num_check = 2; % Число символов Синхронизации
snr = [-26:2:0];

% массивы данных
numinfobits = 10000;
numcodebits = numinfobits;
data = randi([0, 1], 1, numinfobits); % формирование массива бит
num_sym = numcodebits/SF;

% частота, фаза, время
Fmin = 0;              % минимальное значение частоты сигнала
Fmax = 30e6+Fmin;     % максимальное значение частоты сигнала
BW = Fmax-Fmin;        % ширина полосы частот
fc = 2200e6;            % центральная частота
band = BW;
Ts = (2^SF)/BW;        % длительность сигнала


%% ================================= модуляция
os = 1;
fs = BW*os;

[mod_chirp, check_chirp, downch, check_data] = lorax_modified( data, fs, Base, Ts, SF, BW);
N = length(downch);

%% ================================= BER
tic
for n = 1:length(snr)
    n 

    [numErr, NumData] = deal(0);

    while  NumData < 1e5
        
        % АБГШ 
        [rxSig, nvar] = awgn(mod_chirp, snr(n),'measured');

        % демодуляция
        [sv, fourier] = delorax_modified( length(data), SF, downch, rxSig);
        dbits = de2bi(sv,SF).';
        rxbits = dbits(:).';

        % подсчет БЕР с учетом задержки
        err = sum(rxbits~=data);

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
save('snr_lora.mat','snr')


