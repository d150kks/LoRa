clc
clear all
close all

tic
%% ================================= Переменные

% коэффициенты
SF = 9;        % коэффициент расширения спектра (от 7 до 12)
bits2sym = SF-4;
BW = 8.3e6;
BW_actual = 8.3e6;
resamp_factor = BW/BW_actual;
fc = 1000e6;

LORA = myLoRaClass(SF,BW);
Base = LORA.Base;
downch = LORA.downch;
Ts = LORA.Ts;

LORA2 = myLoRaClass(10,BW);
sync_sym = LORA2.downch;

intrlv_state = 13;
num_pre = 4;
warning('off')
% return

% profile on
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
[data_ldpc_intrlvRM, numcodebitsRM, num_symRM, zeros2end, flagRM] = LORA.RM(data_ldpc_intrlv);

%% ================================= CRC coding
[data_crc_ldpc_codeRM] = LORA.codeCRC(data_ldpc_intrlvRM, num_symRM);

%% ================================= Mодуляция
[mod_chirp, check_data, check_data_no_gray] = LORA.lorax_modified( data_crc_ldpc_codeRM, num_symRM, 1);
tx_preamble = repmat(LORA.chirp,1,num_pre);
pre_len = num_pre*Base;
tx_downch = repmat(downch,1,num_pre);
% sync_sym = resample(downch,10,1);
tx_chirp = [sync_sym, tx_downch, tx_preamble, mod_chirp, tx_preamble]; % DOWN, UP, PAYLOAD
tx_length = length(tx_chirp);   


cor = xcorr(tx_chirp, sync_sym);
figure(10)
plot(abs(cor))
%%
return
I = real(tx_chirp);
Q = imag(tx_chirp);
% csvwrite('preamb.csv', sync_sym)
% csvwrite('I.csv', I)
% csvwrite('Q.csv', Q)
save('data.mat','data')
save('preamb.mat','sync_sym')
writematrix( sync_sym.', 'preamb.csv')
writematrix( I.', 'I.csv')
writematrix( Q.', 'Q.csv')




