clc
clear all
close all

%% ================================= Переменные
% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
rc_size = 4;
rc = (SF-rc_size);
BW = 2e6;
snr = [-20:1:0];
nIter = 100;

LORA = myLoRaClass_RSG(SF,BW); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Base = LORA.Base;
downch = LORA.downch;
chirp = LORA.chirp;
Ts = LORA.Ts;

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

if(1)
    num_sym = length(data_ldpc_intrlv)/rc;
else
    num_sym = ceil(length(data_ldpc_intrlv)/SF);
    bits2add = SF-mod(length(data_ldpc_intrlv), SF);
    data_ldpc_intrlv = [data_ldpc_intrlv, zeros(1, bits2add)];
end
% (num_sym+8)*Ts*1000
% numinfobits*(1/((num_sym+8)*Ts))/1024
% return
%% ================================= Mодуляция
% RS
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified_crcrs(data_ldpc_intrlv, num_sym);

% % LORA
% [mod_chirp, check_data, check_no_gray] = LORA.lorax_modified(data_ldpc_intrlv, num_sym, 1);
tx_preamble = repmat(chirp,1,num_pre);
tx_downch = repmat(downch,1,num_pre);
sync_sym = myLoRaClass_RSG(SF+1,BW).downch;

tx_chirp = [sync_sym, tx_downch, tx_preamble, mod_chirp];
tx_length = length(tx_chirp);



save('tx_bits.mat','data')
save('tx_chirp.mat','tx_chirp')
return
%% ================================= BER

% Channel
% h11 = zeros(1, tx_length);
% h11(1) = 1;
% h11(2) = 0.5;
% h11(5) = 0.3;
% h11 = load('h.mat')
h11 = zeros(1, tx_length);
h11(1) = 1;
h11(2) = 0.9;
h11(3) = 0.8;
h11(5) = 0.8;
h11(7) = 0.7;
h11(9) = 0.5;
h11(11) = 0.45;
h11(13) = 0.3;
h11(15) = 0.2;
h11(16) = 0.1;
H11 = fft(h11);

tx_chirp_h = ifft( fft(tx_chirp).*H11 );
% tx_chirp_h = tx_chirp;

% вводим частотный сдвиг
fps = BW/Base;
freq_shift = fps*0.5;
dphi=freq_shift*2*pi*(1/BW);% сдвиг

for j=1:tx_length
    tx_chirp_hf(j)=tx_chirp_h(j)*exp(1i*dphi*j);
end

% % Time delay
% delay = randi(50);
% tx_chirp_hft = [zeros(1,delay), tx_chirp_hf, zeros(1, tx_length)];
tx_chirp_hft = tx_chirp_hf;

tic
% snr = 10;
for n = 1:length(snr)
    fprintf('Iter left: %d\n', length(snr)-n+1) 

    [numErr, NumData] = deal(0);

    for iter = 1:nIter
        
        % АБГШ 
        rxSig = awgn(tx_chirp_hft, snr(n), 'measured');

        % Time sync
%         [rxSig_corr, cor] = LORA.CORRELATION(rxSig, sync_sym, tx_length);
%         rxSig_corr = rxSig_corr(Base*2+1:end);
        rxSig_corr = rxSig(Base*2+1:end);

        % Freq Sync
%         rx_preamble = rxSig_corr(num_pre*Base+1:num_pre*2*Base);
%         corrected_signal = rxSig_corr(num_pre*2*Base+1:end);
        [freq_data, corrected_signal, rx_preamble] = LORA.LORA_FREQ_ESTIM_v3(rxSig_corr, num_pre);

        % Demodulation
        % RS
        [soft_bits, hard_bits, sv_rs, sv, fourier, fourier_rs] = LORA.delorax_crcrs( corrected_signal, num_sym, tx_preamble, rx_preamble);

%         % LORA
%         [soft_bits, hard_bits, sv_decode, sv, fourier] = LORA.delorax_modified( corrected_signal, num_sym, tx_preamble, rx_preamble);
%         hard_bits = hard_bits(1:end-bits2add);
%         soft_bits = soft_bits(1:end-bits2add);


        % LDPC decoding
        maxnumiter = 10;
        rx_data_ldpc = randdeintrlv(soft_bits, intrlv_state);
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
% save('lora_rs_ber.mat','BER')
% save('snr_crc.mat','snr')

