clc
clear all
close all


tic

%% ================================= Переменные

% коэффициенты
SF = 11;        % коэффициент расширения спектра (от 7 до 12)
bits2sym = SF-4;
BW = 30e6;

LORA = myLoRaClass(SF,BW);
Base = LORA.Base;
downch = LORA.downch;

num_pre = 8;

% return
%% ================================= LDPC coding
[cfgLDPCEnc,cfgLDPCDec] = LORA.generateConfigLDPC(3/4, 1944);
        
% Number of message bits
numinfobits = cfgLDPCEnc.NumInformationBits; 
numcodebits = cfgLDPCEnc.BlockLength; 

% Message/Iformation bits
data = randi([0 1],1, numinfobits); 
% data_ldpc_code = data;
data_ldpc_code = ldpcEncode(data.', cfgLDPCEnc).';

%% ================================= Rate matching
[data_ldpc_codeRM, numcodebitsRM, num_symRM, zeros2end, flagRM] = LORA.RM(data_ldpc_code);

% bits2sym = SF-4;
% num_sym = numinfobits/bits2sym;
% bits2end = mod( numinfobits, bits2sym);
% zeros2end = bits2sym-bits2end;
% 
% if( bits2end~=0 )
%     numbitsRM = numinfobits+zeros2end;
%     dataRM = [data, zeros(1,zeros2end)];
%     num_symRM = numbitsRM/bits2sym;
% else
%     numbitsRM = numinfobits;
%     dataRM = data;
%     num_symRM = num_sym;
% end
% return

%% ================================= CRC coding
[data_crc_ldpc_codeRM] = LORA.codeCRC(data_ldpc_codeRM, num_symRM);

%% ================================= Mодуляция
snr = -10;
[mod_chirp, check_data] = LORA.lorax_modified( data_crc_ldpc_codeRM, num_symRM, 1);
tx_preamble = repmat(LORA.chirp,1,num_pre);
tx_sig = [tx_preamble, mod_chirp];
% tx_sig = conv(tx_sig, [0.3, -0.2, 0.6, -0.5], 'same');


%% ================================= АБГШ 
rx_sig = awgn(tx_sig, snr,'measured');
rx_preamble = rx_sig(1:length(tx_preamble));
rxSig = rx_sig(length(tx_preamble)+1:end);



%% ================================= демодуляция
aos = 5;
[soft_bits, hard_bits, sv, sv_cor, fourier] = LORA.DELORAX_CRC( rxSig, num_symRM, tx_preamble, rx_preamble, aos);


%% ================================= декодирование CRC
[data_crc_decodeRM] = LORA.decodeCRC(soft_bits, num_symRM, zeros2end, flagRM);


%% ================================= LDPC decoding
maxnumiter = 100;
rxbits = ldpcDecode(data_crc_decodeRM.', cfgLDPCDec, maxnumiter).';
% rxbits = data_crc_decodeRM;

%% ================================= Debugging
errStats = sum(data~=rxbits);
fprintf('Number of errors = %d\n', errStats)

Tsig = LORA.Ts*num_symRM;
bps = numinfobits/Tsig;

fprintf('T = %.2f\n', Tsig)
fprintf('Bits per Sec = %.2f B/s \n', bps)
fprintf('Bits per Sec = %.2f B/s \n', sum(sv_cor-check_data))
% fprintf('Spectral efficiency = %d\n', errStats)



% figure(1)
% plot(normalize(double(rxbits)))
% hold on
% plot(normalize(data))
% 
% return
