clc
clear all
close all


tic

%% ================================= Переменные
data = load("data.mat").data;

% коэффициенты
SF = 11;        % коэффициент расширения спектра (от 7 до 12)
bits2sym = SF-4;
BW = 30e6;
fc = 2200e6;

LORA = myLoRaClass(SF,BW);
Base = LORA.Base;
downch = LORA.downch;
Ts = LORA.Ts;

num_pre = 4;


%% ================================= LDPC coding
numcodebits = 648;
data_ldpc_code = data;

%% ================================= Rate matching
[data_ldpc_codeRM, numcodebitsRM, num_symRM, zeros2end, flagRM] = LORA.RM(data_ldpc_code);


%% ================================= CRC coding
[data_crc_ldpc_codeRM] = LORA.codeCRC(data_ldpc_codeRM, num_symRM);
% data_crc_ldpc_codeRM = data_ldpc_codeRM;
% num_symRM = numcodebits/SF;


%% ================================= Mодуляция
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified( data_crc_ldpc_codeRM, num_symRM, 1);
tx_preamble = repmat(LORA.chirp,1,num_pre);
tx_chirp = [downch, tx_preamble, mod_chirp]; % DOWN, UP, PAYLOAD
tx_length = length(tx_chirp);
sig_length = tx_length-Base;

 
%% ================================= Correlation

channel_chirp_sto = load("rx_data_in_time.mat").channel_chirp_sto;
[cor,lags] = xcorr(channel_chirp_sto, downch);
cor(round(end/2):end) = cor(round(end/2):end).*gausswin(length(cor(round(end/2):end))).';
[max_amp, max_idx] = max(abs(cor));
start = lags(max_idx);
% start = start-Base*10;
% channel_chirp_corr = channel_chirp_sto(abs(start)+1:abs(start)+tx_length);
channel_chirp_corr = channel_chirp_sto(abs(start)+1:end);


% figure(1)
% plot(real(tx_chirp))
% hold on
% plot(real(channel_chirp_corr))
% N = Base;
% channel_chirp_corr = channel_chirp_corr(Base+1:end);
% for i=1:num_pre
%     fourier_data(i,:) = abs(fft( channel_chirp_corr(i*N-N+1:i*N).*downch ));
% end 

% figure(1)
% plot(real(normalize(tx_chirp(Base+1:end))))
% hold on
% plot(real(normalize(channel_chirp_corr)))
% xlim([1 Base*num_pre])
% 
% figure(2)
% stem(fourier_data.')
% title('fdata')
% 
% [cor2,lags] = xcorr(channel_chirp_sto, tx_preamble);
% figure(3)
% plot(abs(cor2))
% hold on
% plot(abs(cor)*10)
% return




% %% ================================= Frequency correction
% % rx_preamble = channel_chirp_corr(1:Base*num_pre);
% % corrected_signal = channel_chirp_corr(Base*num_pre+1:end);
% 
% % [corrected_signal, rx_preamble] = LORA.CFO(channel_chirp_corr, num_pre);
% [freq_data, corrected_signal, rx_preamble] = LORA.LORA_FREQ_ESTIM(channel_chirp_corr, num_pre, sig_length);
% % [corrected_signal] = LORA_RLS(14, tx_preamble, rx_preamble, corrected_signal);
%% ================================= Frequency correction
N = Base;
fps = BW/Base;

rx_downch = channel_chirp_corr(1:Base);
rx_preamb = channel_chirp_corr(Base+1:Base*(num_pre+1));

% STO AND CFO EST
[r1, r2, fup_idx, fdown_idx, phi1, phi2] = fraq(rx_preamb, rx_downch, downch, Base);
CFOint = fps*(fup_idx-fdown_idx)/2;
CFOfrac = (phi2-phi1)/(2*pi*(N/BW));
STOint = ceil((fup_idx-fdown_idx)/2);
rx_corr_est = channel_chirp_corr((1+Base)-STOint:(0+Base)-STOint+sig_length);
rx_preamb = rx_corr_est(1:Base*num_pre);
rx_signal = rx_corr_est(Base*num_pre+1:end);

channel_chirp_treshold = rx_preamb;

% figure(1)
% stem(abs(r1))
% hold on
% stem(abs(r2))
% title('Coarse')
% fprintf('offset err = %.2f\n', abs(freq_shift-(CFOint+CFOfrac)))
% fprintf('CFOint = %.2f\n', CFOint)
% fprintf('CFOfrac = %.2f\n', CFOfrac)
% fprintf('STOint = %.2f\n\n', STOint)
% return

% 1. Coarse estimation ____________________________________________________
% Определяем Hz/samp
for i = 1:num_pre
    fourier(i,:) = fftshift( abs(fft(  [channel_chirp_treshold(i*N-N+1:N*i).*downch, zeros(1,Base*4)]  )));
    [~, ind1] = max( fourier(i,:) );
    pre_align(i) = (ind1-1)-N*5/2;  
end

% pre_align
% figure(1)
% stem(fourier.')
% return
est1 = (mean(pre_align))*fps/5;
% est1 = round(mean(pre_align))*fps/5;
dphi1 = est1*2*pi*(1/BW); % сдвиг

[channel_chirp_realign] = DPHI_COMP(channel_chirp_treshold, dphi1);
fourier_check_coarse = abs(fft(channel_chirp_realign(1:N).*downch));

% 2. Дробная оценка частоты _______________________________________________
% По двум последовательным одинаковым чирпам как в (1.) вычисляем CFO_fraq
left_ref = downch(1:N/2);
right_ref = downch(N/2+1:N);

for i = 1:num_pre
    chirp2est = channel_chirp_realign(i*N-N+1:i*N);
    left_half = chirp2est(1:N/2);
    right_half = chirp2est(N/2+1:N);
    
    bpf3 = fft( [left_half.*left_ref, zeros(1,N*4/2)]);
    bpf4 = fft( [right_half.*right_ref, zeros(1,N*4/2)]);
    
    [max_a3] = max(bpf3);
    [max_a4] = max(bpf4);
    
    a11=max(max_a3);
    a12=max(max_a4);
    est2_reg(i) = (angle(a12)-angle(a11));
end
est2 = mean(est2_reg)/(pi*Ts);
dphi2 = est2*2*pi*Ts/Base;

[channel_chirp_frac_est] = DPHI_COMP(channel_chirp_realign, dphi2);


% 3. Точное устранение фазового набега ____________________________________
for i = 1:num_pre-1
    argumon(i) = sum(channel_chirp_frac_est(i*N-N+1:N*i).*conj(channel_chirp_frac_est(i*N+1:N*i+N)));
end

[arg] = -mean(angle(argumon));
est3 = arg/(2*pi*Ts);
dphi3 = est3*2*pi/BW; % сдвиг

[channel_chirp_frac_est_v2] = DPHI_COMP(channel_chirp_frac_est, dphi3);
% channel_chirp_frac_est_v2 = channel_chirp_frac_est;
fourier_check_fraq = abs(fft(channel_chirp_frac_est_v2(1:N).*downch));


% 3. Fine estimation ______________________________________________________
for i = 1:num_pre
    fourier2(i,:) = fftshift(abs(fft(  [channel_chirp_frac_est_v2(i*N-N+1:N*i).*downch, zeros(1,Base*4)]  )));
    [~, ind1] = max( fourier2(i,:) );
    pre_align2(i) = (ind1-1)-N*5/2;
end

est4 = (mean(pre_align2))*fps/5;
dphi4 = est4*2*pi*(1/BW); % сдвиг

[cccc] = DPHI_COMP(channel_chirp_frac_est_v2, dphi4);
fourier_check_fine = abs(fft(cccc(1:N).*downch));


figure(1)
subplot(221)
stem(fourier.')
title('None')

subplot(222)
stem(fourier_check_coarse)
title('Coarse')

subplot(223)
stem(fourier_check_fraq)
title('Frac')

subplot(224)
stem(fourier_check_fine)
title('Fine')

fprintf('est1 = %.2f\n', est1)
fprintf('est2 = %.2f\n', est2)
fprintf('est3 = %.2f\n', est3)
fprintf('est4 = %.2f\n\n', est4)

dphi_full = dphi1+dphi2+dphi3+dphi4;
% dphi_full = dphi1+dphi3+dphi4;
% dphi_full = dphi1+dphi3;
% dphi_full = dphi1+dphi2;
% dphi_full = dphi1;
% dphi_full=0;
corrected_signal = zeros(1,length(rx_corr_est));
for j = 1:length(rx_corr_est)
    corrected_signal(j) = rx_corr_est(j).*exp(1i*dphi_full*j*(-1));
end
rx_preamble = corrected_signal(1:num_pre*N);
corrected_signal = corrected_signal(num_pre*N+1:end);


% figure(1)
% subplot(211)
% plot(real(normalize(tx_preamble)))
% hold on
% plot(real(normalize(channel_chirp_realign)))
% legend('ideal', 'err')
% 
% subplot(212)
% plot(abs(cor))
% return
%% ================================= демодуляция
aos = 3;
[soft_bits, hard_bits, sv, sv_cor, fourier] = LORA.DELORAX_CRC( corrected_signal, num_symRM, tx_preamble, rx_preamble, aos);
% [sv_decode, fourier] = LORA.delorax_cyclic(num_symRM, corrected_signal);
% dbits = de2bi(sv_decode, SF)';
% hard_bits = dbits(:)';

%% ================================= декодирование CRC
[data_crc_decodeRM] = LORA.decodeCRC(hard_bits, num_symRM, zeros2end, flagRM);
%         data_crc_decodeRM = hard_bits;


%% ================================= LDPC decoding
%         maxnumiter = 100;
%         rxbits = ldpcDecode(data_crc_decodeRM.', cfgLDPCDec, maxnumiter).';
rxbits = data_crc_decodeRM;




%% ================================= Debugging
bit_err = sum(abs(data~=rxbits));
err_pos = find(check_no_gray~=sv_cor);

%         est1 = freq_data{1};
%         est2 = freq_data{2};
%         est3 = freq_data{3};
%         est_full = freq_data{4};


fprintf('time_shift = %d\n', start)
fprintf('bit_err = %d\n', bit_err)
%         fprintf('est1 = %.2f\n', est1)
%         fprintf('est2 = %.2f\n', est2)
%         fprintf('est3 = %.2f\n', est3)
%         fprintf('est_full = %.2f\n', est_full)

fprintf('check G = %s\n', num2str(check_data(err_pos)))
fprintf('sv G    = %s\n', num2str(sv(err_pos)))
fprintf('check  = %s\n', num2str(check_no_gray(err_pos)))
fprintf('sv_cor = %s\n', num2str(sv_cor(err_pos)))

% fprintf('check  = %s\n', num2str(find(check_no_gray~=sv_cor)))


%% 
tx_chirp2 = tx_chirp(Base*5+1:end);

for i = 1:length(err_pos)
    fft_idl(i,:) = abs(fft(tx_chirp2(err_pos(i)*Base-Base+1: err_pos(i)*Base).*downch));
    fft_err(i,:) = abs(fft(corrected_signal(err_pos(i)*Base-Base+1: err_pos(i)*Base).*downch));

    sig_idl(i,:) = tx_chirp2(err_pos(i)*Base-Base+1: err_pos(i)*Base);
    sig_err(i,:) = corrected_signal(err_pos(i)*Base-Base+1: err_pos(i)*Base);
end

% c = 3;  % big foffset without correction
c= 1;



% figure(1)
% dphi = -0.002;
% for i=1:6
%     % stem( (fft_idl(c,:)))
%     % hold on
% %     dphi = dphi+0.01;
%     for j = 1:length(corrected_signal)
%         corrected_signal(j) = corrected_signal(j).*exp(1i*dphi*j*(-1));
%     end
%     fft_err = abs(fft(corrected_signal(err_pos(1)*Base-Base+1: err_pos(1)*Base).*downch));
% %     stem( (fft_idl(1,:)))
% %     hold on
%     stem( (fft_err))
% %     xlim([170 230])
%     drawnow
%     pause(0.5)
% %     clf
% end
% return


figure(2)
subplot(211)
plot(real(normalize(sig_idl(c,:))))
hold on
plot(real(normalize(sig_err(c,:))))
legend('ideal', 'err')

subplot(212)
stem( (fft_idl(c,:)))
hold on
stem( (fft_err(c,:)))
% xlim([460 510])
% legend('ideal', 'err')


% check = 268  272  312  329  337  341  348  357  359  362  372  375

figure(3)
plot((real(sig_err(c,:))))
hold on
plot(real(channel_chirp_corr))
% legend('hard', 'soft')

function [r1, r2, fup_idx, fdown_idx, phi1, phi2] = fraq(rx_preamb, rx_downch, downch, Base);
    
    r1 = fft( [rx_preamb(1:Base).*downch, zeros(1,Base*4)] );
    r11 = fft(rx_preamb(Base+1:2*Base).*downch);
    r2 = fft( [rx_downch.*conj(downch), zeros(1,Base*4)] );
    [~, fup_idx] = max(abs(r1));
    [~, fdown_idx] = max(abs(r2));
    [~, fup_idx2] = max(abs(r11));
%     if(fup_idx>(Base/2))
%         fup_idx = (fup_idx-1)-Base;
%     end
%###################
    phi2 = angle(r1(fup_idx));
    phi1 = angle(r2(fup_idx2));

    if(fup_idx>(Base*5/2))
        fup_idx = (fup_idx-1)-Base*5;
    end
    if(fup_idx2>(Base*5/2))
        fup_idx2 = (fup_idx2-1)-Base*5;
    end
    if(fdown_idx>(Base*5/2))
        fdown_idx = (fdown_idx-1)-Base*5;
    end

    fup_idx = ceil(fup_idx/5); %###################
    fdown_idx = ceil(fdown_idx/5); %###################
%     fup_idx2 = ceil(fup_idx2/5);
%###################
%     [~, fup_idx2] = max(abs(r11));

% phi2=0;
% phi1=0;

end

        function [output_signal] = DPHI_COMP(input_signal, dphi)
            len = length(input_signal);
            output_signal = zeros(1, len);
            for j=1:len
                output_signal(j)=input_signal(j).*exp(1i*dphi*j*(-1));
            end
        end