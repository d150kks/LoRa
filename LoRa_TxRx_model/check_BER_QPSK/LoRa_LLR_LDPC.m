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
numinfobits = 324;
numcodebits = 648;
data = randi([0, 1], 1, numinfobits); % формирование массива бит
num_sym = numcodebits/SF;

% частота, фаза, время
Fmin = 0;              % минимальное значение частоты сигнала
Fmax = 30e6+Fmin;     % максимальное значение частоты сигнала
BW = Fmax-Fmin;        % ширина полосы частот
fc = 2200e6;            % центральная частота
band = BW;
Ts = (2^SF)/BW;        % длительность сигнала

%% ================================= Mapper
bitmap = de2bi(0:2^SF-1);
M0 = zeros(Base,SF);
M1 = zeros(Base,SF);
A=[];
for nBit=1:SF
    for nSym = 1:Base
        if( bitmap(nSym,nBit) == 0)
            M0(nSym,nBit)= nSym;
        else
            M1(nSym,nBit)= nSym;
        end
    end
end


%% ================================= LDPC coding
P = [ 0 -1 -1 -1  0  0 -1 -1  0 -1 -1  0  1  0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    22  0 -1 -1 17 -1  0  0 12 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1 -1 -1 -1
    6 -1  0 -1 10 -1 -1 -1 24 -1  0 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1 -1 -1
    2 -1 -1  0 20 -1 -1 -1 25  0 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1 -1
    23 -1 -1 -1  3 -1 -1 -1  0 -1  9 11 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1
    24 -1 23  1 17 -1  3 -1 10 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1
    25 -1 -1 -1  8 -1 -1 -1  7 18 -1 -1  0 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1
    13 24 -1 -1  0 -1  8 -1  6 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1 -1 -1
    7 20 -1 16 22 10 -1 -1 23 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1 -1
    11 -1 -1 -1 19 -1 -1 -1 13 -1  3 17 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1
    25 -1  8 -1 23 18 -1 14  9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0
    3 -1 -1 -1 16 -1 -1  2 25  5 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0
    ];

blockSize = 27;
pcmatrix = ldpcQuasiCyclicMatrix(blockSize,P);
cfgLDPCEnc = ldpcEncoderConfig(pcmatrix);
codeword = ldpcEncode(data.', cfgLDPCEnc).';

%% ================================= модуляция
os = 1;
fs = BW*os;

[mod_chirp, check_chirp, downch, check_data] = lorax_modified( codeword, fs, Base, Ts, SF, BW);
N = length(downch);

%% ================================= BER
tic
for n = 1:length(snr)
    n 

    [numErr, NumData] = deal(0);

    while  NumData < 1e4
        
        % АБГШ 
        [rxSig, nvar] = awgn(mod_chirp, snr(n),'measured');

        % демодуляция
        [L, sv, fourier] = delorax_llr( length(codeword), SF, downch, rxSig, M0, M1, nvar);


        % Decoding
        maxnumiter = 50;
        cfgLDPCDec = ldpcDecoderConfig(pcmatrix);
        rxbits = ldpcDecode(L.', cfgLDPCDec, maxnumiter).';

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
save('lora_llr_ldpc_ber.mat','BER')
save('snr_lora.mat','snr')






%%
function [L, sv, fourier] = delorax_llr( length_data, SF, downch, chirp, M0, M1, nvar)

num = length(downch);
% B=2^SF;

num_sym = length_data/SF;
sv = zeros(1,num_sym);
L = [];
for i = 1:num_sym

    d = chirp(num*i-num+1:num*i).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
    
    fourier = abs(fft(d));            % переводим результат в область частот
    [peakMak, indexMax] = max( fourier ); % находим щелчок  частоты в чирпе

    % вычисляем значение кодового слова исходя из базы сигнала
    sv(i) = indexMax-1;

    % LLR

        for nBit=1:SF
          m0 = M0(:,nBit);
          m0(m0==0)=[];
          m1 = M1(:,nBit);
          m1(m1==0)=[];
          LLR = -(1/nvar)*(min( (peakMak-fourier(m0)).^2 ) - min( (peakMak-fourier(m1)).^2 ));
          L = [L LLR];
        end

end

end