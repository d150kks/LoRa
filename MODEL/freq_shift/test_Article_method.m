% _________________________________________________________________________
% On Time-Frequency Synchronization in LoRa System: From
% Analysis to Near-Optimal AlgorithmAnalysis to Near-Optimal Algorithm
% _________________________________________________________________________

clc
clear all
close all

%% ================================= Переменные
% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
rc_size = 4;
rc = (SF-rc_size);
BW = 125e3;
nIter = 1000;

LORA = myLoRaClass_true(SF,BW);
Base = LORA.Base;
downch = LORA.downch;
chirp = LORA.chirp;
ts = LORA.ts;

num_sym = 1;
numinfobits = SF*num_sym;
data = randi([0 1],1, numinfobits); 

%% ================================= Канал (AWGN + Phase shift)
fps = BW/Base;
max_peak_shift = 3;
resamp_factor = 10;

snr = 0;
num_pre_list = 2:2:8;

%% =================================  Modulation
n_pre = 4;

[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified(data, num_sym, 1);
tx_preamble = repmat(chirp,1,n_pre);
tx_downch = repmat(downch,1,n_pre);

tx_chirp = [tx_downch, tx_preamble, mod_chirp];
tx_length = length(tx_chirp);

tx_chirp_fshift = zeros(1,tx_length);
est_err_list = zeros(1,nIter);

tic
for iter = 1:nIter
%     freq_shift = randi([-round(fps*max_peak_shift), round(fps*max_peak_shift)]);
    freq_shift = fps*(-max_peak_shift+(iter*2*max_peak_shift)/nIter);
    dphi = freq_shift*2*pi*ts;% сдвиг
%     dphi = 0;
    
    % вводим частотный сдвиг
    for j=1:tx_length
        tx_chirp_fshift(j)=tx_chirp(j)*exp(1i*dphi*j);
    end


%     time_shift = round(-max_peak_shift + 2*max_peak_shift*rand(), 1);
%     time_shift = round(-max_peak_shift+(iter*2*max_peak_shift)/nIter, 1);
    time_shift = 0;
    tx_chirp_ftshift = resample(tx_chirp_fshift, resamp_factor, 1);
    tx_chirp_ftshift = circshift(tx_chirp_ftshift, time_shift*10);
    tx_chirp_ftshift_dec = resample(tx_chirp_ftshift, 1, resamp_factor);

    rx_preamble = tx_chirp_ftshift_dec(n_pre*Base+1:n_pre*Base*2);
    [freq_data, ~, corrected_preamb] = LORA.LORA_FREQ_ESTIM_v3(tx_chirp_ftshift_dec, n_pre);

    % Algorithm 1
%     [idx, val, fourier] = maxDFT(rx_preamble, downch, 0, n_pre, Base);
%     [index_max, bin_max, fourier_ar] = exhaustiveSynchro(rx_preamble, downch, n_pre, Base);
% index_max

%     figure(1)
% %     plot( fourier )
%     plot( fourier_ar )
%     return

    STO_est = freq_data{1};
    est_full = freq_data{2}+freq_data{3};
    est1 = freq_data{2};
    est2 = freq_data{3};
    STO_reg(iter) = STO_est(1);
    time_shift_reg(iter) = time_shift;


    
    %% =================================  AWGN
    tx_chirp_ftshift_dec_n = awgn(tx_chirp_ftshift_dec, snr, 'measured');

    %% =================================  Freq estimation
    [freq_data, ~, ~] = LORA.LORA_FREQ_ESTIM_v3(tx_chirp_ftshift_dec_n, n_pre);

    STO_est = freq_data{1};
    est_full = freq_data{2}+freq_data{3};
    est1 = freq_data{2};
    est2 = freq_data{3};

    est_err_list(iter) = freq_shift-est1-est2;

end
toc   



        
figure(1); hold on
plot((time_shift_reg), '-')
plot(-(STO_reg))

% figure(2); hold on
% plot(time_shift_reg+STO_reg, '-')
% plot(-())

function [idx, val, fourier] = maxDFT(r, s, k, num_pre, Base)
    fourier = 0;
    r = circshift(r, k);
%     s = circshift(s, k);
    for i=1:num_pre
        fourier = fourier + fftshift(abs(fft(r(Base*i-Base+1:Base*i).*s)));   % перемножаем входной и опорный ОБРАТНый чирп

    end
    [val, idx] = max( fourier );
end

function [index_max, bin_max, fourier_ar] = exhaustiveSynchro(r, s, num_pre, Base)
    index_max = 0;
    bin_max = 0;
    value_max = 0;
    for i=1:Base
        [idx, val, fourier] = maxDFT(r, s, i, num_pre, Base);
        if(val>value_max)
            value_max = val;
            bin_max = idx;
            index_max = i;
        end
    fourier_ar(i,:)=fourier;
    end
end









