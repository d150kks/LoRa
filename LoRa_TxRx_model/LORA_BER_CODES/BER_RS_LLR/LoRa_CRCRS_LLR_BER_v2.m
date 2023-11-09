clc
clear all
close all


%% ================================= Переменные
% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
rc_size = 4;
rc = (SF-rc_size);
BW = 2e6;
snr = [-16:1:0];
nIter = 10;

LORA = myLoRaClass_true(SF,BW);
Base = LORA.Base;
downch = LORA.downch;
chirp = LORA.chirp;

num_pre = 4;
intrlv_state = 13;

% num_sym = 1000;
% numinfobits = num_sym*bits2sym; 
% data = randi([0 1],1, numinfobits); 


%% ================================= LDPC coding
[cfgLDPCEnc,cfgLDPCDec] = LORA.generateConfigLDPC(5/6, 1944);
        
% Number of message bits
numinfobits = cfgLDPCEnc.NumInformationBits; 
numcodebits = cfgLDPCEnc.BlockLength; 

% Message/Iformation bits
data = randi([0 1],1, numinfobits); 
data_ldpc = ldpcEncode(data.', cfgLDPCEnc).';
data_ldpc_intrlv = randintrlv(data_ldpc, intrlv_state);
num_sym = length(data_ldpc_intrlv)/rc;

%% ================================= Mодуляция
% [mod_chirp, check_data, check_no_gray] = LORA.lorax_modified_crcrs(data, num_sym);
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified_crcrs(data_ldpc_intrlv, num_sym);
tx_preamble = repmat(chirp,1,num_pre);
tx_downch = repmat(downch,1,num_pre);
sync_sym = myLoRaClass_true(SF+1,BW).downch;

%% ================================= BER

% Channel
h11 = zeros(1, length(mod_chirp));
h11(1) = 1;
h11(2) = 0.5;
h11(5) = 0.3;
H11 = fft(h11);

mod_chirp = ifft( fft(mod_chirp).*H11 );
% mod_chirp = mod_chirp;

% вводим частотный сдвиг
fps = BW/Base;
freq_shift = fps*5.5;
dphi=freq_shift*2*pi*(1/BW);% сдвиг

for j=1:length(mod_chirp)
    mod_chirp(j)=mod_chirp(j)*exp(1i*dphi*j);
end

% Time delay
delay = randi(50);
tx_chirp = [sync_sym, tx_downch, tx_preamble, mod_chirp];
tx_length = length(tx_chirp);
tx_chirp = [zeros(1,delay), tx_chirp, zeros(1, Base)];

tic
snr = -8;
for n = 1:length(snr)
    fprintf('Iter: %d\n', n) 

    [numErr, NumData] = deal(0);

    for iter = 1:nIter
        
        % АБГШ 
        rxSig = awgn(tx_chirp, snr(n), 'measured');
        [rxSig_corr, cor] = LORA.CORRELATION(rxSig, sync_sym, tx_length);
        rxSig_corr = rxSig_corr(Base*2+1:end);

        [freq_data, corrected_signal, rx_preamble] = LORA.LORA_FREQ_ESTIM_v3(rxSig_corr, num_pre);
        [rx_preamble, var] = awgn(tx_preamble,snr(n),'measured');
%         var
%         noiseVar = 10.^(-snr(n)/10)
%         nvar = std( abs(normalize(tx_preamble)-((rx_preamble-mean(tx_preamble))./std(tx_preamble))).^2 )

        % демодуляция
%         [soft_bits, hard_bits, sv_rs, sv, fourier, fourier_rs] = LORA.delorax_crcrs( rxSig, num_sym);
        [soft_bits, hard_bits, sv_rs, sv, fourier, fourier_rs] = LORA.delorax_crcrs( rxSig, num_sym, tx_preamble, rx_preamble);

        % LDPC decoding
        maxnumiter = 100;
        rx_data_ldpc = randdeintrlv((soft_bits), intrlv_state);
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
% save('lora_crcrs_ber.mat','BER')
% save('lora_rs_ber2.mat','BER')
% save('snr_crc.mat','snr')

