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
LORA.OS = 1;

num_sym = 1;
numinfobits = SF*num_sym;
data = randi([0 1],1, numinfobits); 

%% ================================= Канал (AWGN + Phase shift)
fps = BW/Base;
max_peak_shift = 3;

snr = [-12: 0];
num_pre_list = 2:8;

tic
snr = 5;
for n = 1:length(snr)
    fprintf('Iter left: %d\n', length(snr)-n+1)

    %% =================================  Modulation
    n_pre = 8;
    [mod_chirp, check_data, check_no_gray] = LORA.lorax_modified(data, num_sym, 1);
    tx_preamble = repmat(chirp,1,n_pre);
    
    tx_chirp = [tx_preamble, mod_chirp];
    tx_length = length(tx_chirp);

    tx_chirp_fshift = zeros(1,tx_length);
    est_err_list = zeros(1,nIter);

    for iter = 1:nIter
        freq_shift = randi([-round(fps*max_peak_shift), round(fps*max_peak_shift)]);
        freq_shift = fps*0.5;
        dphi = freq_shift*2*pi*ts;% сдвиг
        
        % вводим частотный сдвиг
        for j=1:tx_length
            tx_chirp_fshift(j)=tx_chirp(j)*exp(1i*dphi*j);
        end
        
        %% =================================  AWGN
        tx_chirp_fshift_n = awgn(tx_chirp_fshift, snr(n), 'measured');

        %% =================================  Freq estimation
%         [freq_data, rx_preamb_comp3] = LORA.LORA_FREQ_ESTIM_ghan(tx_chirp_fshift_n, n_pre);
        [freq_data, rx_preamb_comp3] = LORA.LORA_FREQ_ESTIM_old(tx_chirp_fshift_n, n_pre);
    
        freq_data
        est_full = freq_data{1};
        est1 = freq_data{2};
        est2 = freq_data{3};
        est3 = freq_data{4};
        est_err_list(iter) = freq_shift-est1-est2-est3;
        freq_shift-est1-est2-est3
return
    end
    est_err_from_pre(n) = std(est_err_list);
end

toc   

figure(1)
histogram(est_err_list, 60)

% figure(2)
% stem(est_err_list)
% ylabel('Frequency, Hz');
% grid on

figure(3)
plot(snr, est_err_from_pre)
ylabel('Frequency, Hz');
grid on


% save('fest_gss_snr.mat','est_err_from_pre'); % сохранение неразвернутой преамбулы
% save('fest_ghanaatian_snr.mat','est_err_from_pre'); % сохранение неразвернутой преамбулы
save('fest_our_snr.mat','est_err_from_pre'); % сохранение неразвернутой преамбулы



