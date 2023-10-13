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


%% ================================= LDPC coding
maxnumiter = 100;
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
        data_code_ldpc = ldpcEncode(data.', cfgLDPCEnc).';
        % data_code_ldpc = data;

        
        %% ================================= модуляция
        mod_sig = qammod(data_code_ldpc.', M, 'InputType','bit').';
        
        
        %% ================================= BER
        tic
        for n = 1:length(snr)

            [numErr, NumData] = deal(0);
        
            while  NumData < 1e5
                
                % АБГШ 
                rxSig = awgn(mod_sig,snr(n),'measured');
        
                % демодуляция
                demod_sig = qamdemod(rxSig.', M, 'OutputType','bit').';
        
                % декодирование LDPC
                data_decode_ldpc = double(~ldpcDecode((2*demod_sig-1).', cfgLDPCDec, maxnumiter).');
                
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
save('qpsk_config_ldpc_ber.mat','BER')
save('snr_config_ldpc.mat','snr')

