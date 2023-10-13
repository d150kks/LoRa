clc
clear all
close all

%% ================================= Переменные

% коэффициенты
SF = 8;        % коэффициент расширения спектра (от 7 до 12)
Base = 2^SF;   % База сигнала
num_sym = 2000;% Число передаваемых символов
num_pre = 8;   % Число символов Преамбулы
num_check = 2; % Число символов Синхронизации
snr = [-20:2:0];

% массивы данных
numbits = SF*num_sym;          % число бит
data = randi([0, 1], 1, numbits); % формирование массива бит

% частота, фаза, время
Fmin = 0;              % минимальное значение частоты сигнала
Fmax = 30e6+Fmin;     % максимальное значение частоты сигнала
BW = Fmax-Fmin;        % ширина полосы частот
fc = 1000e6;            % центральная частота
band = 30e6;
Ts = (2^SF)/BW;        % длительность сигнала


%% ================================= модуляция
os = 1;
fs = BW*os;
ts=1/fs;

[mod_chirp, check_chirp, downch, check_data] = lorax_modified( data, fs, Base, Ts, SF, BW);
N = length(downch);

% Forming Signal
tx_chirp = mod_chirp;
tx_length = length(tx_chirp);


%% ================================= BER
tic
for n = 1:length(snr)

    [numErr, NumData] = deal(0);

    while  NumData < 1e5
        % АБГШ демодуляция и декодирование
        rxSig = awgn(mod_chirp,snr(n),'measured');
        [sv, fourier] = delorax_modified( length(data), SF, downch, rxSig);
        dbits = de2bi(sv,SF)';
        dbits = dbits(:)';
        
        
        % подсчет БЕР с учетом задержки
        err = sum(dbits~=data);

        % Increment the error and bit counters
        numErr = numErr + err;        
        NumData = NumData + numbits;
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

% save('lora_ber.mat','BER')
% save('lora_snr.mat','snr')
