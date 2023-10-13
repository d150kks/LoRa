clc
clear all
close all


tic

%% ================================= Переменные

% коэффициенты
SF = 9;        % коэффициент расширения спектра (от 7 до 12)
Base = 2^SF;   % База сигнала
num_sym = 100;% Число передаваемых символов
num_pre = 8;   % Число символов Преамбулы
num_check = 2; % Число символов Синхронизации
snr = [-36:2:-16];

% массивы данных
numinfobits = 2310;          % число бит
data = randi([0, 1], 1, numinfobits); % формирование массива бит
crc_type = 4;
bits_per_sym = (SF-crc_type);
num_sym = numinfobits/bits_per_sym;

% частота, фаза, время
Fmin = 0;              % минимальное значение частоты сигнала
Fmax = 30e6+Fmin;     % максимальное значение частоты сигнала
BW = Fmax-Fmin;        % ширина полосы частот
fc = 2200e6;            % центральная частота
band = BW;
Ts = (2^SF)/BW;        % длительность сигнала

%% ================================= CRC coding
data_win = zeros(1, bits_per_sym);
for i=1:num_sym
    data_win = data(i*bits_per_sym-bits_per_sym+1:i*bits_per_sym);
    CRC_Bits = CRC4(data_win.').';
    data_code(i*SF-SF+1:i*SF) = [data_win, CRC_Bits].';
end


%% ================================= модуляция
os = 1;
fs = BW*os;
ts=1/fs;

% ##################
data_resh = reshape(data_code, SF, []);
data_rep = repmat(data_resh, 4, 1);
data_rep = reshape(data_rep, 1, []);
% ##################

[mod_chirp, check_chirp, downch, check_data] = lorax_modified( data_rep, fs, Base, Ts, SF, BW);
N = length(downch);


%% ================================= BER
tic
for n = 1:length(snr)
    n 

    [numErr, NumData] = deal(0);

    while  NumData < 1e5
        
        % АБГШ 
        rxSig = awgn(mod_chirp,snr(n),'measured');

        % демодуляция
        [sv, fourier] = delorax_modified2( length(data_rep), SF, downch, rxSig);
        dbits = de2bi(sv,SF).';
        dbits = dbits(:).';

        % декодирование CRC
        data_decode = zeros(1, numinfobits);
        for i=1:num_sym
            data_win = dbits(i*SF-SF+1:i*SF);
            data_decode(i*bits_per_sym-bits_per_sym+1:i*bits_per_sym) = data_win(bits_per_sym-bits_per_sym+1:bits_per_sym);
        end

%         figure(1)
%         plot(data-data_decode)
%         return
        % подсчет БЕР с учетом задержки
        err = sum(data_decode~=data);

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
save('lora_crc_dd_ber.mat','BER')
save('snr_crc.mat','snr')

%%

function [sv, fourier] = delorax_modified2( length_data, SF, downch, chirp)
num = length(downch);
B=2^SF;
% css = [];
for i = 1:(length_data/SF)/4

pair_chirp = chirp((i-1)*num*4+1: num*4*(i));
d = pair_chirp(1:num).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
d2 = pair_chirp(num+1:num*2).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
d3 = pair_chirp(num*2+1:num*3).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
d4 = pair_chirp(num*3+1:num*4).*downch;   % перемножаем входной и опорный ОБРАТНый чирп

fourier = abs(fft(d));            % переводим результат в область частот
fourier2 = abs(fft(d2));            % переводим результат в область частот
fourier3 = abs(fft(d3));            % переводим результат в область частот
fourier4 = abs(fft(d4));            % переводим результат в область частот
[~, indexMax] = max( (fourier+fourier2+fourier3+fourier4) ); % находим щелчок  частоты в чирпе
% css = fourier
%%%%%%%
% вычисляем значение кодового слова исходя из базы сигнала

sv(i) = indexMax-1;
% indexMax-1;
% if indexMax>B
%     sv(i) = B - (num - indexMax) - 1;
% else
%     sv(i) = indexMax - 1;
% end

end

end


