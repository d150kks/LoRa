clc
clear all
close all


tic

%% ================================= Переменные

% коэффициенты
SF = 12;        % коэффициент расширения спектра (от 7 до 12)
bits2sym = SF-4;
BW = 125e3;
% fc = 2220e6;

LORA = myLoRaClass(SF,BW);
Base = LORA.Base;
downch = LORA.downch;
Ts = LORA.Ts;

intrlv_state = 13;
num_pre = 4;
warning('off')

% profile on
%% ================================= LDPC coding
[cfgLDPCEnc,cfgLDPCDec] = LORA.generateConfigLDPC(1/2, 1944);
        
% Number of message bits
numinfobits = cfgLDPCEnc.NumInformationBits; 
numcodebits = cfgLDPCEnc.BlockLength; 

% Message/Iformation bits
data = randi([0 1],1, numinfobits); 
% data = zeros(1,numinfobits);
data_ldpc = ldpcEncode(data.', cfgLDPCEnc).';
data_ldpc_intrlv = randintrlv(data_ldpc, intrlv_state);

%% ================================= Rate matching
[data_ldpc_intrlvRM, numcodebitsRM, num_symRM, zeros2end, flagRM] = LORA.RM(data_ldpc_intrlv);

%% ================================= CRC coding
[data_crc_ldpc_codeRM] = LORA.codeCRC(data_ldpc_intrlvRM, num_symRM);

%% ================================= Mодуляция
[mod_chirp, check_data, check_data_no_gray] = LORA.lorax_modified( data_crc_ldpc_codeRM, num_symRM, 1);
tx_preamble = repmat(LORA.chirp,1,num_pre);
pre_len = num_pre*Base;
tx_downch = repmat(downch,1,num_pre);
sync_sym = resample(downch,2,1);
tx_chirp = [sync_sym, tx_downch, tx_preamble, mod_chirp, tx_preamble]; % DOWN, UP, PAYLOAD
tx_length = length(tx_chirp);

% return
% channel_chirp_sto = rx_data_in_time.';
channel_chirp_sto = [zeros(1, 13134),tx_chirp];

%% ================================= Correlation
[channel_chirp_corr, cor] = LORA.CORRELATION(channel_chirp_sto, sync_sym, tx_length, 1);
channel_chirp_corr = channel_chirp_corr(Base*2+1:end);

% input_signal = channel_chirp_sto;
% if( log2(tx_length)<ceil(log2(tx_length)) )
%     padding = 2^ceil(log2(tx_length))-tx_length;
% else
%     padding=0;
% end
% 
% input_signal_fft = fft([input_signal, zeros(1,padding)]);
% sync_sym_fft = fft(sync_sym, length(input_signal_fft));
% cor = ifft(input_signal_fft.*conj(sync_sym_fft));
% [~, start] = max(abs(cor));
% channel_chirp_corr = channel_chirp_sto(start:start+tx_length-1);

figure(1); hold on;
% plot(real(channel_chirp_sto))
% plot(real(channel_chirp_corr))
plot(abs(cor))

% length(channel_chirp_sto)+padding
% 524288
return
%% ================================= Frequency correction
[freq_data, corrected_signal, rx_preamble] = LORA.LORA_FREQ_ESTIM_v3(channel_chirp_corr, num_pre);
[corrected_signal] = LORA.STO_COMP(corrected_signal, num_pre);

%% ================================= демодуляция
[soft_bits, hard_bits, sv, sv_cor, fourier] = LORA.DELORAX_CRC( corrected_signal, num_symRM, tx_preamble, rx_preamble);


%% ================================= декодирование CRC
[rx_data_ldpc_intrlv] = LORA.decodeCRC(soft_bits, num_symRM, zeros2end, flagRM);


%% ================================= LDPC decoding
maxnumiter = 100;
rx_data_ldpc = randdeintrlv(rx_data_ldpc_intrlv, intrlv_state);
rx_data = ldpcDecode(rx_data_ldpc.', cfgLDPCDec, maxnumiter).';

%% ================================= Debugging
bit_err = sum(abs(data~=rx_data));

STO = freq_data{1};
est1 = freq_data{2};
est2 = freq_data{3};
freq_shift_est = est1+est2;

fprintf('Num Err: %d\n', bit_err)
fprintf('STO:     %.2f\n\n', STO)

fprintf('Offset:    %.2f\n', entered_offset)
fprintf('Est 1:    %.2f\n', est1)
fprintf('Est 2:    %.2f\n', est2)
fprintf('Offset Est:  %.2f\n', freq_shift_est)
fprintf('Offset Err:  %.2f\n\n', abs(freq_shift_est+entered_offset))










