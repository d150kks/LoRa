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

num_sym = 1000;
numinfobits = SF*num_sym;
data = randi([0 1],1, numinfobits); 
num_pre = 4;

%% =================================  Modulation
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified(data, num_sym, 1);
tx_preamble = repmat(chirp,1,num_pre);
tx_downch = repmat(downch,1,num_pre);

tx_chirp = [tx_downch, tx_preamble, mod_chirp];
tx_length = length(tx_chirp);


%% ================================= Канал (AWGN + Phase shift)
fps = BW/Base;
max_peak_shift = 3;
snr = [-16:0];

tic
for n = 1:length(snr)

    fprintf('Iter left: %d\n', length(snr)-n) 

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
        
        %% =================================  AWGN
        tx_chirp_fshift_n = awgn(tx_chirp_fshift, snr(n), 'measured');

        %% =================================  Freq estimation
        [freq_data, corrected_signal, rx_preamb_comp3] = LORA.LORA_FREQ_ESTIM_v3(tx_chirp_fshift_n, num_pre);
    
        STO_est = freq_data{1};
        est_full = freq_data{2}+freq_data{3};
        est1 = freq_data{2};
        est2 = freq_data{3};
        est_err_list(iter) = freq_shift-est1-est2;

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

figure(3)
semilogy(snr, BER)
% ylabel('Frequency, Hz');
grid on



% save('plots/fest_ber_new_std.mat','BER'); % сохранение неразвернутой преамбулы



