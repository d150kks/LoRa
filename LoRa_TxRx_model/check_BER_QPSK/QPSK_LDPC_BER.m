clc
clear all
close all


tic

%% ================================= Переменные

% коэффициенты
M = 4;
snr = [-20:2:20];

numinfobits = 486;          % число бит
data = randi([0, 1], 1, numinfobits); % формирование массива бит

%% ================================= LDPC coding
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
codeword = ldpcEncode(data.', cfgLDPCEnc).';

%% ================================= модуляция
mod_sig = qammod(codeword.', M, 'InputType','bit').';

%% ================================= BER
tic
for n = 1:length(snr)
    n 

    [numErr, NumData] = deal(0);

    while  NumData < 1e5
        
        % АБГШ 
        rxSig = awgn(mod_sig,snr(n),'measured');

        % демодуляция
%         demod_sig = qamdemod(rxSig.', M, OutputType='approxllr').';
        demod_sig = qamdemod(rxSig.', M, 'OutputType','bit').';

%         figure(1)
%         plot(demod_sig)
%         hold on
%         plot(demod_sig2)
%         return

        % Decoding
        maxnumiter = 50;
        cfgLDPCDec = ldpcDecoderConfig(pcmatrix);
%         rxbits = ldpcDecode(demod_sig.', cfgLDPCDec, maxnumiter).';
        rxbits = ldpcDecode(100*(1-2*demod_sig).', cfgLDPCDec, maxnumiter).';

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
save('qpsk_ldpc_ber.mat','BER')
% save('snr.mat','snr')

