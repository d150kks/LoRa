clc
clear all
close all


tic
% rng default
% rng(13)
%% ================================= Переменные

% коэффициенты
SF = 9;        % коэффициент расширения спектра (от 7 до 12)
bits2sym = SF-4;
BW = 125e3;
fc = 2200e6;

LORA = myLoRaClass(SF,BW);
Base = LORA.Base;
Ts = LORA.Ts;
downch = LORA.downch;

num_pre = 8;
numcodebits = 648;
data = randi([0 1],1, numcodebits); 
data_ldpc_code = data;

%% ================================= Rate matching
[data_ldpc_codeRM, numcodebitsRM, num_symRM, zeros2end, flagRM] = LORA.RM(data_ldpc_code);


%% ================================= CRC coding
[data_crc_ldpc_codeRM] = LORA.codeCRC(data_ldpc_codeRM, num_symRM);


%% ================================= Mодуляция
[mod_chirp, check_data] = LORA.lorax_modified( data_crc_ldpc_codeRM, num_symRM, 1);
tx_preamble = repmat(LORA.chirp,1,num_pre);

tx_chirp = [downch, tx_preamble, mod_chirp];
tx_length = length(tx_chirp);


%% ================================= АБГШ 
snr = -10;
% freq_shift = 400;
freq_shift = 244*18.5;
ts = (1/BW);
dphi=freq_shift*2*pi*ts;% сдвиг

% вводим частотный сдвиг
for j=1:tx_length
    shift_sig(j)=tx_chirp(j)*exp(1i*dphi*j);
end
% shift_sig = tx_chirp;

% Time Delay
delay = randi([10,100]);
% delay = 0;
rx_sig = awgn( [zeros(1,delay), shift_sig, zeros(1,1e3)], snr, 'measured');
% rx_sig = awgn( shift_sig, snr, 'measured');


%% ================================= Correlation
[cor,lags] = xcorr(rx_sig, downch);
[max_amp, max_idx] = max(abs(cor));
start = lags(max_idx);
% start = delay;


rx_corr = rx_sig(abs(start)+1:abs(start)+tx_length);
rx_downch = rx_corr(1:Base);
rx_preamb = rx_corr(Base+1:Base*(num_pre+1));
rx_signal = rx_corr(Base+1:end);
% rx_signal = rx_corr(Base*(num_pre+1)+1:end);

% figure(1)
% plot(abs(cor))
%% ================================= 1. Preamble Synchronization
N = Base;
fps = BW/Base;



%% ================================= Frequency correction

[freq_data_1, output_signal_est_1, rx_preamb_est_1] = LORA.CFO(rx_corr, num_pre);
[freq_data_2, output_signal_est_2, rx_preamb_est_2] = LORA.LORA_FREQ_ESTIM(rx_corr, num_pre);
STOint_1 = freq_data_1{1};
CFOint_1 = freq_data_1{2};
% CFOfraq_1 = freq_data_1{3};
CFOfraq_1 = 0;
FEraw_1 = freq_data_1{4};

STOint_2 = freq_data_2{1};
CFOint_2 = freq_data_2{2};
CFOfraq_2 = freq_data_2{3};
FEraw_2 = freq_data_2{4};

% Debugging
fprintf('delay     = %d\n', delay)
fprintf('delay est = %d\n', start)
start_cor_1 = start-STOint_1;
start_cor_2 = start-STOint_2;
fprintf('delay_err = %d; %d\n\n', abs(delay-start_cor_1), abs(delay-start_cor_2))

fprintf('STOint = %d; %d\n', STOint_1, STOint_2)
fprintf('CFOint  = %.2f, %.2f\n', CFOint_1, CFOint_2)
fprintf('CFOfrac = %.2f, %.2f\n', CFOfraq_1, CFOfraq_2)
fprintf('FEraw   = %.2f, %.2f\n\n', FEraw_1, FEraw_2)

fprintf('freq offset = %.2f\n', freq_shift)
fer1 = abs(freq_shift-(CFOint_1+CFOfraq_1+FEraw_1));
fer2 = abs(freq_shift-(CFOint_2+CFOfraq_2+FEraw_2));
fprintf('freq offset err = %.2f; %.2f\n', fer1, fer2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


return
%% ================================= демодуляция
[sv_decode, sv, fourier] = LORA.delorax_modified(input_signal_est_full, num_symRM);
dbits = de2bi(sv_decode, SF)';
hard_bits = dbits(:)';
% aos = 5;
% [soft_bits, hard_bits, sv, sv_cor, fourier] = LORA.DELORAX_CRC( rx_payload, num_symRM, tx_preamble, rx_preamble, aos);

err_pos = find(check_data~=sv);
fprintf('check  = %s\n', num2str(check_data(err_pos)))
fprintf('sv = %s\n', num2str(sv(err_pos)))

figure(1)
stem(abs(fourier),'LineWidth',2)
grid on
%% ================================= декодирование CRC
[data_crc_decodeRM] = LORA.decodeCRC(hard_bits, num_symRM, zeros2end, flagRM);
%         data_crc_decodeRM = hard_bits;


%% ================================= LDPC decoding
%         maxnumiter = 100;
%         rxbits = ldpcDecode(data_crc_decodeRM.', cfgLDPCDec, maxnumiter).';
rxbits = data_crc_decodeRM;




%% ================================= Debugging
bit_err = sum(abs(data~=rxbits));
% err_pos = find(check_no_gray~=sv_cor);

%         est1 = freq_data{1};
%         est2 = freq_data{2};
%         est3 = freq_data{3};
%         est_full = freq_data{4};

fprintf('bit_err = %d\n', bit_err)



function [r1, r2, fup_idx, fdown_idx, phi1, phi2] = fraq(rx_preamb, rx_downch, downch, Base);
    
    r1 = fft(rx_preamb(1:Base).*downch);
    r11 = fft(rx_preamb(Base+1:2*Base).*downch);
    r2 = fft(rx_downch.*conj(downch));
    [~, fup_idx] = max(abs(r1));
    [~, fdown_idx] = max(abs(r2));
    if(fdown_idx>(Base/2))
        r2(fdown_idx)=0;
        [~, fdown_idx] = max(abs(r2));
    end

    [~, fup_idx2] = max(abs(r11));
    phi2 = angle(r1(fup_idx));
    phi1 = angle(r2(fup_idx2));

end