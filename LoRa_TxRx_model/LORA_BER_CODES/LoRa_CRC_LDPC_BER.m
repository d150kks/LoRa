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
snr = [-20:2:-0];

% массивы данных
numinfobits = 486;
numcodebits = 648;
data = randi([0, 1], 1, numinfobits); % формирование массива бит
% data = zeros(1,numinfobits);
crc_type = 4;
bits_per_sym = (SF-crc_type);
num_sym = numcodebits/bits_per_sym;


% частота, фаза, время
Fmin = 0;              % минимальное значение частоты сигнала
Fmax = 30e6+Fmin;     % максимальное значение частоты сигнала
BW = Fmax-Fmin;        % ширина полосы частот
fc = 2200e6;            % центральная частота
band = 30e6;
Ts = (2^SF)/BW;        % длительность сигнала

%% ================================= LDPC coding
maxnumiter = 100;
P = [16 17 22 24  9  3 14 -1  4  2  7 -1 26 -1  2 -1 21 -1  1  0 -1 -1 -1 -1
     25 12 12  3  3 26  6 21 -1 15 22 -1 15 -1  4 -1 -1 16 -1  0  0 -1 -1 -1
     25 18 26 16 22 23  9 -1  0 -1  4 -1  4 -1  8 23 11 -1 -1 -1  0  0 -1 -1
      9  7  0  1 17 -1 -1  7  3 -1  3 23 -1 16 -1 -1 21 -1  0 -1 -1  0  0 -1
     24  5 26  7  1 -1 -1 15 24 15 -1  8 -1 13 -1 13 -1 11 -1 -1 -1 -1  0  0
      2  2 19 14 24  1 15 19 -1 21 -1  2 -1 24 -1  3 -1  2  1 -1 -1 -1 -1  0
    ];
blockSize = 27;
pcmatrix = ldpcQuasiCyclicMatrix(blockSize,P);

cfgLDPCEnc = ldpcEncoderConfig(pcmatrix);
cfgLDPCDec = ldpcDecoderConfig(pcmatrix);

data_code_ldpc = ldpcEncode(data.', cfgLDPCEnc).';
% data_code_ldpc = data;

%% ================================= CRC coding
data_win = zeros(1, bits_per_sym);
for i=1:num_sym
    data_win = data_code_ldpc(i*bits_per_sym-bits_per_sym+1:i*bits_per_sym);
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

    while  NumData < 1e5
        
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

        % декодирование LDPC
        data_decode_ldpc = double(~ldpcDecode((2*data_decode-1).', cfgLDPCDec, maxnumiter).');
        
        % подсчет БЕР с учетом задержки
        err = sum(data_decode_ldpc~=data);

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
% 
save('lora_crc_ldpc_ber.mat','BER')
save('snr.mat','snr')

