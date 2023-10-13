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
% tx_sync = repmat(LORA.sync,1,2);
tx_chirp = [downch, tx_preamble, mod_chirp];
% tx_chirp = [downch, tx_preamble, tx_sync, mod_chirp]; %%%% New Preamble
tx_length = length(tx_chirp);
sig_length = tx_length-Base;


%% ================================= АБГШ 
snr = -8;
% freq_shift = -244*32.5;
freq_shift = 244*11.5;
% freq_shift = -50;
ts = (1/BW);
dphi=freq_shift*2*pi*ts;% сдвиг

% вводим частотный сдвиг
for j=1:tx_length
    shift_sig(j)=tx_chirp(j)*exp(1i*dphi*j);
end
% shift_sig = tx_chirp;
coeffs = [-0.000576496332561738	0.0356137014750986	0.0348073675130581	-0.0999613697894719	0.0360018880205920	0.526437991958931	0.526437991958931	0.0360018880205920	-0.0999613697894719	0.0348073675130581	0.0356137014750986	-0.000576496332561738];
shift_sig = conv(shift_sig, coeffs, 'same');

% figure(1)
% plot(real(shift_sig))
% return

% Time Delay
delay = randi([100,1000]);
% delay = 0;
rx_sig = awgn( [zeros(1,delay), shift_sig, zeros(1,2e4)], snr, 'measured');
% rx_sig = awgn( shift_sig, snr, 'measured');


%% ================================= Correlation
[cor,lags] = xcorr(rx_sig, downch);
[max_amp, max_idx] = max(abs(cor));
start = lags(max_idx);
% start = delay;


% rx_corr = rx_sig(abs(start)+1:abs(start)+tx_length);
rx_corr = rx_sig(abs(start)+1:end);
rx_downch = rx_corr(1:Base);
rx_preamb = rx_corr(Base+1:Base*(num_pre+1));
% rx_signal = rx_corr(Base+1:end);
% rx_signal = rx_corr(Base*(num_pre+1)+1:end);

% figure(1)
% plot(abs(cor))
% return
%% ================================= 1. Preamble Synchronization
N = Base;
fps = BW/Base;

% spec= 1000*ones(1,100);
% time = fft(spec);
% sig = conv(tx_preamble, time, 'same')/90;
% sig = filter(time, 1,tx_preamble)/90;
% 
% figure(1)
% plot(real(sig))
% hold on
% plot(real(tx_preamble))
% return

%% ================================= Frequency correction
% STO AND CFO EST
[r1, r2, fup_idx, fdown_idx, phi1, phi2] = fraq(rx_preamb, rx_downch, downch, Base);
CFOint = fps*(fup_idx-fdown_idx)/2;
CFOfrac = (phi2-phi1)/(2*pi*(N/BW));
STOint = ceil((fup_idx-fdown_idx)/2);

rx_corr_est = rx_corr((1+Base)-STOint:(0+Base)-STOint+sig_length);
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


channel_chirp_realign = zeros(1,length(channel_chirp_treshold));
for j=1:length(channel_chirp_treshold)
    channel_chirp_realign(j)=channel_chirp_treshold(j).*exp(1i*dphi1*j*(-1));
end
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


channel_chirp_frac_est = zeros(1,length(channel_chirp_realign));
for j=1:length(channel_chirp_realign)
    channel_chirp_frac_est(j) = channel_chirp_realign(j).*exp(1i*dphi2*j*(-1));
end

% figure(1)
% stem(abs(bpf3))
% est2
% return

% 3. Точное устранение фазового набега ____________________________________
for i = 1:num_pre-1
    argumon(i) = sum(channel_chirp_frac_est(i*N-N+1:N*i).*conj(channel_chirp_frac_est(i*N+1:N*i+N)));
end

[arg] = -mean(angle(argumon));
est3 = arg/(2*pi*Ts);
dphi3 = est2*2*pi/BW; % сдвиг

% точное устранение фазового набега
for j=1:length(channel_chirp_frac_est)
    channel_chirp_frac_est_v2(j) = channel_chirp_frac_est(j).*exp(1i*dphi3*j*(-1));
end
fourier_check_fraq = abs(fft(channel_chirp_frac_est_v2(1:N).*downch));


% return
% 3. Fine estimation ______________________________________________________
for i = 1:num_pre
    fourier2(i,:) = fftshift(abs(fft(  [channel_chirp_frac_est(i*N-N+1:N*i).*downch, zeros(1,Base*4)]  )));
    [~, ind1] = max( fourier2(i,:) );
    pre_align2(i) = (ind1-1)-N*5/2;
end

% est3 = round(mean(pre_align2))*fps/5;
est4 = (mean(pre_align2))*fps/5;
dphi4 = est3*2*pi*(1/BW); % сдвиг

cccc = zeros(1,length(channel_chirp_frac_est));
for j=1:length(channel_chirp_frac_est)
    cccc(j)=channel_chirp_frac_est(j).*exp(1i*dphi4*j*(-1));
end

% Debugging _______________________________________________________________
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

fprintf('delay     = %d\n', delay)
fprintf('delay est = %d\n', start)
start_cor = start-STOint;
fprintf('delay_err = %d\n\n', abs(delay-start_cor))

fprintf('offset = %.2f\n', freq_shift+0)
fprintf('offset err = %.2f\n', abs(freq_shift-(est1+est2+est3)))
fprintf('est1 = %.2f\n', est1)
fprintf('est2 = %.2f\n', est2)
fprintf('est3 = %.2f\n\n', est3)

% return
% dphi_full = dphi1+dphi2+0+0;
dphi_full = dphi1+dphi2+dphi3+0;
input_signal_est_full = zeros(1,length(rx_corr_est));
for j = 1:length(rx_corr_est)
    input_signal_est_full(j) = rx_corr_est(j).*exp(1i*dphi_full*j*(-1));
end
input_signal_est_full = input_signal_est_full(num_pre*N+1:end);

% pre_align
% figure(1)
% % plot(pre_align)
% stem(abs(fourier))
% return
% 
% % точное устранение фазового набега
% channel_chirp_frac_est = zeros(1,length(channel_chirp_treshold));
% for j=1:length(channel_chirp_realign)
%     channel_chirp_frac_est(j) = channel_chirp_realign(j).*exp(1i*dphi2*j*(-1));
% end
% 
% % 3. fine estimation
% % Устранение фазового сдвига
% for i = 1:num_pre-1
%     argumon(i*N-N+1:N*i) = channel_chirp_frac_est(i*N-N+1:N*i).*conj(channel_chirp_frac_est(i*N+1:N*i+N));
% end
% 
% [arg] = max(sum(argumon));
% est3 = -angle(arg)/(2*pi*obj.Ts);
% 
% % Устраняем CFO_fraq сигнала
% dphi3 = est3*2*pi*obj.Ts/obj.Base; % сдвиг
% 
% % точное устранение фазового набега
% preamb_align = zeros(1,length(channel_chirp_treshold));
% for j=1:length(channel_chirp_frac_est)
%     preamb_align(j) = channel_chirp_frac_est(j).*exp(1i*dphi3*j*(-1));
% end


% Debugging
% fprintf('delay     = %d\n', delay)
% fprintf('delay est = %d\n', start)
% start_cor = start-STOint;
% fprintf('delay_err = %d\n\n', abs(delay-start_cor))

% fprintf('STOint = %d\n', STOint)
% fprintf('CFOint  = %.2f, %.3f\n', CFOint, CFOint_dphi)
% fprintf('CFOfrac = %.2f, %.3f\n', CFOfrac, CFOfrac_dphi)
% fprintf('FEraw   = %.2f, %.3f\n\n', FEraw, DEraw)
% 
% fprintf('offset = %.2f\n', freq_shift)
% fprintf('offset err = %.2f\n', abs(freq_shift-(CFOint+CFOfrac+FEraw)))


% [freq_data, input_signal_est_full, rx_preamble] = LORA.LORA_FREQ_ESTIM(rx_signal, num_pre);
% est1 = freq_data{1};
% est2 = freq_data{2};
% est3 = freq_data{3};
% est_full = freq_data{4};
% 
% fprintf('time_shift = %d\n', start)
% fprintf('est1 = %.2f\n', est1)
% fprintf('est2 = %.2f\n', est2)
% fprintf('est3 = %.2f\n', est3)
% fprintf('est_full = %.2f\n', est_full)

%% ================================= демодуляция
% [sv_decode, sv, fourier] = LORA.delorax_modified(input_signal_est_full, num_symRM);
[sv_decode, fourier] = LORA.delorax_cyclic(num_symRM, input_signal_est_full);
dbits = de2bi(sv_decode, SF)';
hard_bits = dbits(:)';
% aos = 5;
% [soft_bits, hard_bits, sv, sv_cor, fourier] = LORA.DELORAX_CRC( input_signal_est_full, num_symRM, tx_preamble, rx_preamble, aos);

% err_pos = find(check_data~=sv);
% fprintf('check  = %s\n', num2str(check_data(err_pos)))
% fprintf('sv = %s\n', num2str(sv(err_pos)))

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




% %% ================================= 3. Golden rectangle
% acc = 0.01;     % Интервал оценки
% k=1:Base;    % Массив для ускорения расчетов
% g = 0.618;   % Число золотого сечения
% fps = fps*2;
% for n=1:num_pre
%     % выбираем 
%     r = channel_chirp_realign(n*Base-Base+1:Base*n); % один чирп преамбулы
%     L=-0.5*fps;
%     R=0.5*fps;
%     
%     while (R-L)>acc
%         df1 = L+(1-g)*(R-L);
%         df2 = L+g*(R-L);
%          
%          % формируем 8 опорных чирпов содержахи космпенсирующую частоту
%         Local_chirp1 = exp(-1i*2*pi*(df1*(k+Base))/BW);
%         Local_chirp2 = exp(-1i*2*pi*(df2*(k+Base))/BW);
%         
%         % Компенсируем у преамбулы сдвиг, находим функцию Z(f) в точке df
%         Zi1 = r.*Local_chirp1;
%         Zi2 = r.*Local_chirp2;
%         
%         % Находим функцию F(f) в точках df
%         Ri1 =  abs(sum((Zi1.*downch)));
%         Ri2 =  abs(sum((Zi2.*downch)));
%         
%         if Ri1<Ri2
%             L=df1;
%         end
%         if Ri1>Ri2
%             R=df2;
%         end
%     end
%     est_sqr(n) = (R+L)/2;
% end
% % est_sqr
% est2 = sum(est_sqr)/num_pre;
% dphi2 = est2*2*pi*Ts/Base;

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
