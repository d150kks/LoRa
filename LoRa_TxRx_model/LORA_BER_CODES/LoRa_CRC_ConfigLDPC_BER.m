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
snr = [-20:2:0];

% частота, фаза, время
Fmin = 0;              % минимальное значение частоты сигнала
Fmax = 30e6+Fmin;     % максимальное значение частоты сигнала
BW = Fmax-Fmin;        % ширина полосы частот
fc = 2200e6;            % центральная частота
band = 30e6;
Ts = (2^SF)/BW;        % длительность сигнала

%% ================================= LDPC coding
maxnumiter = 50;
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
        crc_type = 4;
        bits_per_sym = (SF-crc_type);
        num_sym = numcodebits/bits_per_sym;
        
        
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
        N = length(downch);
        
        
        %% ================================= BER
        tic
        for n = 1:length(snr)

            [numErr, NumData] = deal(0);
        
            while  NumData < 1e4
                
                % АБГШ 
                [rxSig, nvar] = awgn(mod_chirp,snr(n),'measured');
        
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
                data_decode_ldpc = double(~ldpcDecode((1/nvar)*(2*data_decode-1).', cfgLDPCDec, maxnumiter).');
                
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
save('lora_crc_config_ldpc_ber.mat','BER')
save('snr_config_ldpc.mat','snr')

