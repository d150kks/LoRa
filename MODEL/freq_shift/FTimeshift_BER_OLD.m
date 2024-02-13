clc
clear all
close all

%% ================================= Переменные
% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
rc_size = 4;
rc = (SF-rc_size);
BW = 125e3;
nIter = 100;

LORA = myLoRaClass_true(SF,BW);
Base = LORA.Base;
downch = LORA.downch;
chirp = LORA.chirp;
ts = LORA.ts;
LORA.OS = 1;

num_sym = 2000;
numinfobits = SF*num_sym;
data = randi([0 1],1, numinfobits); 
num_pre = 8;

%% =================================  Modulation
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified(data, num_sym, 1);
tx_preamble = repmat(chirp,1,num_pre);

tx_chirp = [tx_preamble, mod_chirp];
tx_length = length(tx_chirp);


%% ================================= Канал (AWGN + Phase shift)
fps = BW/Base;
max_peak_shift = 3;
resamp_factor = 30;
snr = [-16:0];

tic
for n = 1:length(snr)

    fprintf('Iter left: %d\n', length(snr)-n+1) 

    [numErr, NumData] = deal(0);
    tx_chirp_fshift = zeros(1,tx_length);
    est_err_list = zeros(1,nIter);

    for iter = 1:nIter
        freq_shift = randi([-round(fps*max_peak_shift), round(fps*max_peak_shift)]);
    %     freq_shift = normrnd(0, fps*max_peak_shift);
        dphi = freq_shift*2*pi*ts;% сдвиг
        
        % вводим частотный сдвиг
        for j=1:tx_length
            tx_chirp_fshift(j)=tx_chirp(j)*exp(1i*dphi*j);
        end

        time_shift = round(-max_peak_shift + 2*max_peak_shift*rand(), 1);
        tx_chirp_ftshift = resample(tx_chirp_fshift, resamp_factor, 1);
        tx_chirp_ftshift = circshift(tx_chirp_ftshift, time_shift*10);
        tx_chirp_ftshift_dec = resample(tx_chirp_ftshift, 1, resamp_factor);
        
        %% =================================  AWGN
        tx_chirp_fshift_n = awgn(tx_chirp_ftshift_dec, snr(n), 'measured');

        %% =================================  Freq estimation
        [freq_data, rx_preamb_comp3, corrected_signal] = LORA.LORA_FREQ_ESTIM_old(tx_chirp_fshift_n, num_pre);
    
        est_full = freq_data{1};
        est1 = freq_data{2};
        est2 = freq_data{3};
        est3 = freq_data{4};

        est_err_list(iter) = freq_shift-est1-est2-est3;

        %% =================================  Demodulation
        % демодуляция
        [soft_bits, hard_bits, sv_decode, sv, fourier] = LORA.delorax_modified( corrected_signal, num_sym, tx_preamble, rx_preamb_comp3);

        % подсчет БЕР с учетом задержки
        err = sum(hard_bits~=data);

        % Increment the error and bit counters
        numErr = numErr + err;        
        NumData = NumData + numinfobits;
    end

    % Estimate the BER for both methods
    BER(n) = numErr/NumData;
    est_err_from_pre(n) = std(est_err_list);

end

toc   

figure(1)
histogram(est_err_list, 60)

figure(2)
stem(est_err_list)
ylabel('Frequency, Hz');
grid on

figure(3)
semilogy(snr, BER)
% ylabel('Frequency, Hz');
grid on



% save('plots/ftest_ber_old_std.mat','BER'); % сохранение неразвернутой преамбулы


