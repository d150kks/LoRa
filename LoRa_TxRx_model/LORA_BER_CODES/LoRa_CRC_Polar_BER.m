clc
clear all
close all


tic

%% ================================= Переменные

% коэффициенты
SF = 8;        % коэффициент расширения спектра (от 7 до 12)
Base = 2^SF;   % База сигнала
% num_sym = 400;% Число передаваемых символов
num_pre = 8;   % Число символов Преамбулы
num_check = 2; % Число символов Синхронизации
snr = [-30:2:-10];

% массивы данных
numinfobits = 384;
numcodebits = 512;
data = randi([0, 1], 1, numinfobits); % формирование массива бит
% data = zeros(1,numinfobits);
crc_type = 4;
bits_per_sym = (SF-crc_type);
num_sym = numcodebits/bits_per_sym;
% numinfobits = SF*num_sym;          % число бит


% частота, фаза, время
Fmin = 0;              % минимальное значение частоты сигнала
Fmax = 30e6+Fmin;     % максимальное значение частоты сигнала
BW = Fmax-Fmin;        % ширина полосы частот
fc = 2200e6;            % центральная частота
band = 30e6;
Ts = (2^SF)/BW;        % длительность сигнала
% return
%% ================================= Polar coding
E = 512;
data_code_polar = nrPolarEncode(data.',E,10,false).';


%% ================================= CRC coding
data_win = zeros(1, bits_per_sym);
for i=1:num_sym
    data_win = data_code_polar(i*bits_per_sym-bits_per_sym+1:i*bits_per_sym);
    CRC_Bits = CRC4(data_win.').';
    data_code(i*SF-SF+1:i*SF) = [data_win, CRC_Bits].';
end

%% ================================= модуляция
os = 1;
fs = BW*os;
ts=1/fs;

[mod_chirp, check_chirp, downch, check_data] = lorax_modified( data_code, fs, Base, Ts, SF, BW);
% [mod_chirp, check_chirp, downch, check_data] = LORAX_CRC( data_code, fs, Base, Ts, SF, BW);
N = length(downch);

% Forming Signal
preamble = repmat(check_chirp,1,num_pre);
%%tx_chirp = [preamble, mod_chirp];
tx_chirp = [downch, preamble, mod_chirp];
tx_length = length(tx_chirp);

%% ================================= BER
tic
for n = 1:length(snr)
    n 

    [numErr, NumData] = deal(0);

    while  NumData < 1e4
        
        % АБГШ 
        rxSig = awgn(mod_chirp,snr(n),'measured');

        % демодуляция
        aos = 5;
        [demod_bits, sv_cor, sv, fourier] = DELORAX_CRC( length(data_code), SF, downch, rxSig, aos);

        % декодирование CRC
        data_decode = zeros(1, numcodebits);
        for i=1:num_sym
            data_win = demod_bits(i*SF-SF+1:i*SF);
            data_decode(i*bits_per_sym-bits_per_sym+1:i*bits_per_sym) = data_win(bits_per_sym-bits_per_sym+1:bits_per_sym);
        end

        % декодирование Polar
        L = 8;
        data_decode_polar = nrPolarDecode((2*~data_decode-1).',numinfobits,E,L,10,false,24).';
        
        % подсчет БЕР с учетом задержки
        err = sum(data_decode_polar~=data);

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
% save('lora_crc_polar_ber.mat','BER')
% save('snr.mat','snr')