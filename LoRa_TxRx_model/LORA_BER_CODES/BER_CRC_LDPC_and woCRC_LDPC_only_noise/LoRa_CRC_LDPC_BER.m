clc
clear all
close all


tic

%% ================================= Переменные
% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
bits2sym = SF-4;
BW = 30e6;
snr = [-16:1:0];
nIter = 500;

LORA = myLoRaClass(SF,BW);
Base = LORA.Base;
downch = LORA.downch;

num_pre = 8;

intrlv_state = 13;
% nbits = 486; 
% data = randi([0 1],1, nbits); 

%% ================================= LDPC coding
[cfgLDPCEnc,cfgLDPCDec] = LORA.generateConfigLDPC(3/4, 648);
        
% Number of message bits
numinfobits = cfgLDPCEnc.NumInformationBits; 
numcodebits = cfgLDPCEnc.BlockLength; 

% Message/Iformation bits
data = randi([0 1],1, numinfobits); 
data_ldpc = ldpcEncode(data.', cfgLDPCEnc).';
data_ldpc_intrlv = randintrlv(data_ldpc, intrlv_state);

%% ================================= Rate matching
[data_ldpc_intrlvRM, numcodebitsRM, num_sym, zeros2end, flagRM] = LORA.RM(data_ldpc_intrlv);

%% ================================= CRC coding
[data_ldpc_intrlvRM_crc] = LORA.codeCRC(data_ldpc_intrlvRM, num_sym);

% return
%% ================================= Mодуляция
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified( data_ldpc_intrlvRM_crc, num_sym, 1);

num_pre = 8;
tx_preamble = repmat(LORA.chirp,1,num_pre);
tx_chirp = [tx_preamble, mod_chirp]; % DOWN, UP, PAYLOAD

%% ================================= BER

tic

for n = 1:length(snr)
    n 

    [numErr, NumData] = deal(0);

    for iter = 1:nIter
        
        % АБГШ 
        rxSig = awgn(tx_chirp,snr(n),'measured');
        rx_preamble = rxSig(1:Base*num_pre);
        rx_sig = rxSig(Base*num_pre+1:end);

        % демодуляция
        aos = 3;
        [soft_bits, hard_bits, sv, sv_cor, fourier] = LORA.DELORAX_CRC(rx_sig, num_sym, tx_preamble, rx_preamble, aos);

        % CRC decoding
        [rx_data_ldpc_intrlv] = LORA.decodeCRC(soft_bits, num_sym, zeros2end, flagRM);

        % LDPC decoding
        maxnumiter = 100;
        rx_data_ldpc = randdeintrlv(rx_data_ldpc_intrlv, intrlv_state);
        rx_data = ldpcDecode(rx_data_ldpc.', cfgLDPCDec, maxnumiter).';

        % подсчет БЕР с учетом задержки
        err = sum(rx_data~=data);

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
save('lora_crc_ldpc_ber_648.mat','BER')
% save('lora_crc_ldpc_ber_1944.mat','BER')
% save('snr_crc.mat','snr')

