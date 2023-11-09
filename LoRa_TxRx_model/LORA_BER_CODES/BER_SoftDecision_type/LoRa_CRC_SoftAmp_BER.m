clc
clear all
close all

% return
tic

%% ================================= Переменные
% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
bits2sym = SF-4;
BW = 30e6;
snr = [-16:1:0];
nIter = 20;

LORA = myLoRaClass_test(SF,BW);
Base = LORA.Base;
downch = LORA.downch;

num_pre = 8;
intrlv_state = 13;

% num_sym = 10;
% nbits = 1200; 
% data = randi([0 1],1, nbits); 
% data = [0 0 1];


%% ================================= LDPC coding
[cfgLDPCEnc,cfgLDPCDec] = LORA.generateConfigLDPC(3/4, 1944);
        
% Number of message bits
numinfobits = cfgLDPCEnc.NumInformationBits; 
numcodebits = cfgLDPCEnc.BlockLength; 

% Message/Iformation bits
data = randi([0 1],1, numinfobits); 
% data = zeros(1,numinfobits);
data_ldpc = ldpcEncode(data.', cfgLDPCEnc).';
data_ldpc_intrlv = randintrlv(data_ldpc, intrlv_state);


%% ================================= Rate matching
[dataRM, numcodebitsRM, num_sym, zeros2end, flagRM] = LORA.RM(data_ldpc_intrlv);


%% ================================= CRC coding
[data_code] = LORA.codeCRC(dataRM, num_sym);

%% ================================= Mодуляция
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified( data_code, num_sym, 1);
tx_preamble = repmat(LORA.chirp,1,num_pre);

%% ================================= BER

tic

h11 = zeros(1, length(mod_chirp));
h11(1) = 1;
h11(2) = 0.5;
h11(5) = 0.3;
H11 = fft(h11);

mod_chirp = ifft( fft(mod_chirp).*H11 );

fps = BW/Base;
freq_shift = fps*0.5; %%%%%%%%%%%%%%%%%%%%%%%%
dphi=freq_shift*2*pi*(1/BW);% сдвиг

% вводим частотный сдвиг
for j=1:length(mod_chirp)
    mod_chirp(j)=mod_chirp(j)*exp(1i*dphi*j);
end


% snr = 10;
for n = 1:length(snr)
    n 

    [numErr, NumData] = deal(0);

    for iter = 1:nIter
        
        % АБГШ 
        rxSig = awgn(mod_chirp,snr(n),'measured');
        rx_preamble = awgn(tx_preamble, snr(n), 'measured');

        % демодуляция
        [soft_bits, hard_bits, sv, sv_cor, fourier] = LORA.DELORAX_CRC( rxSig, num_sym, tx_preamble, rx_preamble);

        % декодирование CRC
        [rx_data_ldpc_intrlv] = LORA.decodeCRC(soft_bits, num_sym, zeros2end, flagRM);
        
        % LDPC decoding
        maxnumiter = 10;
        rx_data_ldpc = randdeintrlv(rx_data_ldpc_intrlv, intrlv_state);
        data_decode = ldpcDecode(rx_data_ldpc.', cfgLDPCDec, maxnumiter).';

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
% save('lora_crc_softAMP_ber.mat','BER')
% save('snr_crc.mat','snr')
figure(2); hold on
plot( -soft_bits)
plot( (2*hard_bits-1).*max(soft_bits))
