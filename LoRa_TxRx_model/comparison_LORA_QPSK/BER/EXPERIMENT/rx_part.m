clc
clear all
close all
return
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

intrlv_state = 13;
num_pre = 4;
warning('off')
% return

num_symRM = 389;
zeros2end = 1;
flagRM = 1;
tx_length = 207360;

[cfgLDPCEnc,cfgLDPCDec] = LORA.generateConfigLDPC(3/4, 1944);
tx_preamble = repmat(LORA.chirp,1,num_pre);
pre_len = num_pre*Base;
LORA2 = myLoRaClass(SF+2,BW);
sync_sym = LORA2.downch;
load('data.mat', '-mat');

%% ================================= RX
load_data = load('rx_data.mat', '-mat');
rx_data = complex(double(load_data.I), double(load_data.Q));
% rx_data = ifft(fftshift(fft(rx_data)));
% rx_data = resample(rx_data,1,2);
% 
% channel_chirp_sto = rx_data_in_time.';
% %         channel_chirp_sto = resample(channel_chirp_sto,1,resamp_factor);

%% ================================= Correlation
[channel_chirp_corr, cor] = LORA.CORRELATION(rx_data, sync_sym, tx_length, 1);
channel_chirp_corr = channel_chirp_corr(Base*4+1:end);


% figure(1)
% plot(imag(rx_data))
% figure(10)
% plot(abs(cor))
% return
%% ================================= Frequency correction
[freq_data, corrected_signal, rx_preamble] = LORA.LORA_FREQ_ESTIM_v3(channel_chirp_corr, num_pre);
%         [corrected_signal, w] = MY_NLMS(Base/2, 0.1, tx_preamble, rx_preamble, corrected_signal);
[corrected_signal] = LORA.STO_COMP(corrected_signal, num_pre);

N = Base;
rx_pre = channel_chirp_corr(pre_len+1:pre_len*2);
for i=1:num_pre
    fourier_rx(i,:) = abs(fft( rx_pre(i*N-N+1:i*N).*downch ));
end
for i=1:num_pre
    fourier_fcor(i,:) = abs(fft( rx_preamble(i*N-N+1:i*N).*downch ));
end
for i=1:num_symRM
    fourier_data(i*N-N+1:i*N) = abs(fft( corrected_signal(i*N-N+1:i*N).*downch ));
end

figure(1)
subplot(221)
stem(fourier_rx.')
title('rx')

subplot(222)
stem(fourier_fcor.')
title('fcor')

subplot(223)
stem(fourier_data.')
title('fdata')

subplot(224)
semilogy( (abs(fftshift(fft(channel_chirp_corr)))).')
ylim([1 1000])
title('spec')
drawnow

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

% fprintf('Offset:    %.2f\n', entered_offset)
fprintf('Est 1:    %.2f\n', est1)
fprintf('Est 2:    %.2f\n', est2)
fprintf('Offset Est:  %.2f\n', freq_shift_est)
fprintf('Offset Err:  %.2f\n\n', abs(freq_shift_est+0))
fprintf('BER:     %.5f\n\n', bit_err/1458) 
