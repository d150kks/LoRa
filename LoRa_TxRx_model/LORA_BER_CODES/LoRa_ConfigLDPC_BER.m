clc
clear all
close all


tic

%% ================================= Переменные
% Class initialization
LORA = myLoRaClass;

% коэффициенты
SF = 8;        % коэффициент расширения спектра (от 7 до 12)
Base = 2^SF;   % База сигнала
% num_sym = 400;% Число передаваемых символов
num_pre = 8;   % Число символов Преамбулы
num_check = 2; % Число символов Синхронизации
snr = [-20:1:0];

% частота, фаза, время
Fmin = 0;              % минимальное значение частоты сигнала
Fmax = 30e6+Fmin;     % максимальное значение частоты сигнала
BW = Fmax-Fmin;        % ширина полосы частот
fc = 2200e6;            % центральная частота
band = 30e6;
Ts = (2^SF)/BW;        % длительность сигнала

%% ================================= LDPC coding
maxnumiter = 10;
rate = [1/2 2/3 3/4 5/6];
cdwlen = [648, 1296, 1944]; % Codeword length
for nLen = 1:3
    nLen
    for nRate = 1:4
        [cfgLDPCEnc,cfgLDPCDec] = generateConfigLDPC(rate(nRate), cdwlen(nLen));
        
        % Number of message bits
        numinfobits = cfgLDPCEnc.NumInformationBits; 
        numcodebits = cfgLDPCEnc.BlockLength; 
        
        % Message/Iformation bits
        data = randi([0 1],1, numinfobits); 
        
        % массивы данных
        num_sym = numcodebits/SF;
        
%         num_sym
        data_code_ldpc = ldpcEncode(data.', cfgLDPCEnc).'; 
        
        %% ================================= модуляция
        os = 1;
        fs = BW*os;
        
        [mod_chirp, check_chirp, downch, check_data] = LORA.lorax_modified( data_code_ldpc, fs, Base, Ts, SF, BW);
        tx_preamble = repmat(check_chirp,1,num_pre);
        tx_sig = [tx_preamble, mod_chirp];
        N = length(downch);

        %% ================================= BER
        tic
        for n = 1:length(snr)

            [numErr, NumData] = deal(0);
        
            while  NumData < 2e4
                
                % АБГШ 
                rx_sig = awgn(tx_sig,snr(n),'measured');
                rx_preamble = rx_sig(1:length(tx_preamble));
                rxSig = rx_sig(length(tx_preamble)+1:end);
        
                % демодуляция
                [L, sv, fourier] = LORA.delorax_llr(num_sym, SF, downch, rxSig, tx_preamble, rx_preamble);
        
                % декодирование LDPC
                data_decode_ldpc = ldpcDecode(L.', cfgLDPCDec, maxnumiter).';
                
                % подсчет БЕР с учетом задержки
                err = sum(data_decode_ldpc~=data);
        
                % Increment the error and bit counters
                numErr = numErr + err;        
                NumData = NumData + numinfobits;
            end
        
            % Estimate the BER for both methods
            BER(n, nLen, nRate) = numErr/NumData;
        end
        toc
    end
end

%%
figure(1)
semilogy(snr,BER(:,3,3),'-*','color','k');
grid
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
hold off
title('SNR');
% 
save('lora_config_ldpc_ber.mat','BER')
save('snr.mat','snr')

