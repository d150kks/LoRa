clc
clear all
close all


tic
% rng default
% rng(13)
%% ================================= Переменные

% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
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
data = zeros(1, numcodebits);
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
sig_length = tx_length-Base;


%% ================================= АБГШ 

% ~~~~~~~~ Freq multipath ~~~~~~~~
h11 = zeros(1, Base);
h11(1) = 1;
h11(2) = 0.5;
h11(5) = 0.3;
H11 = fft(h11);
% H11 = ones(1, Base);

for i=1:tx_length/Base
    tx_sig_fft = fft(tx_chirp(i*Base-Base+1:Base*i));
    tx_chirp2(i*Base-Base+1:Base*i) = ifft(tx_sig_fft.*H11);
end

% ~~~~~~~~ Freq ~~~~~~~~
fps = BW/Base;
snr = 20;
freq_shift = -fps*2.5;
% freq_shift = fps*0.5;
% freq_shift = 0;
ts = (1/BW);
dphi=freq_shift*2*pi*ts;% сдвиг

% вводим частотный сдвиг
for j=1:tx_length
    shift_sig(j)=tx_chirp2(j)*exp(1i*dphi*j);
end
% shift_sig = tx_chirp;

% ~~~~~~~~ Time Delay ~~~~~~~~ 
delay = randi([100,1000]);
delay = 0;
% delay = 90;
pq = 100;
shift_sig = resample(shift_sig,pq,1);
rx_sig = awgn( [zeros(1,delay), shift_sig, zeros(1,2e4)], snr, 'measured');
% rx_sig = resample(rx_sig,1,pq-1);
% rx_sig = rx_sig(48:pq:end);
rx_sig = rx_sig(1:pq:end);
% rx_sig = awgn( shift_sig, snr, 'measured');


%% ================================= Correlation
[cor,lags] = xcorr(rx_sig, downch);
[max_amp, max_idx] = max(abs(cor));
start = lags(max_idx);
start = 3;

rx_corr = rx_sig(abs(start)+1:end);
rx_downch = rx_corr(1:Base);
rx_preamb = rx_corr(Base+1:Base*(num_pre+1));


%% ================================= Timing Correction
% Check current offset
product_arr = [];
for num = 9:100
    [sv, product] = LORA.delorax_cyclic(1, rx_corr(Base*num+1:Base*(num+1)));
    product_arr = [product_arr, product];

    [~, product_max] = max(product);
    if(product_max>Base/2)
        product_max_arr(num) = product_max-Base-1;
    else
        product_max_arr(num) = product_max-1;
    end
end


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sto estimation ~~~~~~~~~~~~~~~~~~~~~~~~~~~
N = Base;
OS = 10;
r1=0;
for i = 1:num_pre
    r1 = r1 + abs( fft([rx_preamb(i*N-N+1:N*i).*downch, zeros(1,N*(OS-1))]) );
end

r2 = fft( [rx_downch.*conj(downch), zeros(1,Base*(OS-1))] );
[~, fup_idx] = max(abs(r1));
[~, fdown_idx] = max(abs(r2));

if(fup_idx>(Base*OS/2))
    fup_idx = (fup_idx-1)-Base*OS;
end
if(fdown_idx>(Base*OS/2))
    fdown_idx = (fdown_idx-1)-Base*OS;
end

fup_idx = (fup_idx/OS); 
fdown_idx = (fdown_idx/OS); 

STO = ((fup_idx-fdown_idx)/2);
if(STO>0)
    STOint = floor(STO);
else
    STOint = ceil(STO);
end
STOfraq = STO-STOint;
fprintf(['STO:     %.2f\n' ...
         'STOint:  %.2f\n' ...
         'STOfraq: %.2f\n'], STO, STOint, STOfraq)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~ STO Compensation ~~~~~~~~~~~~~~~~~~~~~~~~~~~
sto_comp_vec = ones(1,Base)*STOfraq;
if(STO>0)
    sto_comp_vec(Base-STOint+1:end) = 0;
else
    sto_comp_vec(abs(STOint)+1:end) = 0;
end

rx_corr_payload = rx_corr(Base+1:end);
% rx_corr_payload = tx_chirp(Base+1:end);
% rx_corr_payload = rx_corr;
counter = 1;

% for j=1:length(rx_corr_payload)
for j=1:length(tx_chirp)-Base
%     rxSig_sto_comp(j) = rx_corr_payload(j).*exp(1i*2*pi * (-STOfraq/N) );
%     rxSig_sto_comp(j) = rx_corr_payload(j).*exp(1i*2*pi * (-STOint*j/N) );
%     rxSig_sto_comp(j) = rx_corr_payload(j).*exp(1i*2*pi * (-STOint*counter/N) );
%     rxSig_sto_comp(j) = rx_corr_payload(j).*exp(1i*2*pi * (-STOint*j/N) ).*exp(1i*2*pi* (sto_comp_vec(counter)*j/N));
%     rxSig_sto_comp(j) = rx_corr_payload(j).*exp(1i*2*pi * (-STOint*j/N) ).*exp(1i*2*pi * (-STOfraq*j/N) );
    rxSig_sto_comp(j) = rx_corr_payload(j).*exp(1i*2*pi * (-STO*j/N) );
%     check_sig(j) = rx_corr_payload(j)./tx_chirp(Base+j);
%     check_sig2(j) = exp(1i*2*pi * (STOint*j/N) ).*exp(1i*2*pi * (-STOfraq/N) );
%     check_sig2(j) = exp(1i*2*pi * (STOint*j/N +0.19) ).*exp(1i*2*pi* 0.19 );
    counter=counter+1;
    if(counter>Base)
        counter=1;
    end
end
% rxSig_sto_comp = rx_corr_payload;
[freq_data, rxSig_sto_comp, rx_preamble] = LORA.LORA_FREQ_ESTIM_v2(rxSig_sto_comp, num_pre, sig_length/2);

% plot(real(check_sig))
% hold on
% plot(real(check_sig2))
% xlim([1 Base*4])
% return

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~ STO Check ~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Check current offset
product_arr2 = [];
for num = 1:100
    [sv, product2] = LORA.delorax_cyclic(1, rxSig_sto_comp(Base*num+1:Base*(num+1)));
    product_arr2 = [product_arr2, product2];
end

            % ~~~~~~~~ 3. fine estimation ~~~~~~~~ 
            % Устранение фазового сдвига
            N = Base;
            argumon = zeros(1, num_pre-1);
            for i = 1:num_pre-1
                argumon(i) = sum(rxSig_sto_comp(i*N-N+1:N*i).*conj(rxSig_sto_comp(i*N+1:N*i+N)));
            end
            
            [arg] = -mean(angle(argumon));
            est3 = arg/(2*pi*Ts);
            dphi3 = est3*2*pi/BW; % сдвиг

[est0, dphi0] = LORA.COARSE_FREQ_ESTIM(rxSig_sto_comp, num_pre);

fprintf('Est 1:  %.2f\n', est0)
fprintf('Est 3:  %.2f\n', est3)


figure(1)
plot(real(tx_chirp(Base+1:end)))
hold on
plot(real(rxSig_sto_comp))
xlim([1 Base*4])
% return

figure(4)
stem(product_arr)
hold on
stem(product_arr2)
xlim([1 Base*8])

% figure(4)
% subplot(211)
% stem(product_arr)
% subplot(212)
% stem(product_max_arr)
return


























%% ================================= 1. Preamble Synchronization
N = Base;
fps = BW/Base;


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

% fourier_check = abs(fft(rx_preamb(1:N).*downch));
% figure(1)
% stem(abs(r1))
% return
% ~~~~~~~~ TE estimation ~~~~~~~~
% Find TE
% for mp = 1:num_pre
%     fourier = abs(fft(channel_chirp_treshold(mp*Base-Base+1:Base*mp).*downch));
%     [TEraw] = LORA.TIME_ERR(fourier);
%     TEraw_array(mp) = TEraw;
% end
% TEraw_mean = mean( abs(TEraw_array) );
% 
% FEraw = -((BW*TEraw_mean)/Base)/2.685;
% DEraw = FEraw*2*pi*(1/BW);
% 
% % Compensate in preamb
% comp_sig_1 = zeros(1,length(channel_chirp_treshold));
% comp_sig_2 = zeros(1,length(channel_chirp_treshold));
% for j=1:length(channel_chirp_treshold)
%     comp_sig_1(j) = channel_chirp_treshold(j).*exp(1i*DEraw*j*(-1));
%     comp_sig_2(j) = channel_chirp_treshold(j).*exp(1i*DEraw*j);
% end
% 
% for mp = 1:num_pre
%     fourier1 = abs(fft(comp_sig_1(mp*Base-Base+1:Base*mp).*downch));
%     fourier2 = abs(fft(comp_sig_2(mp*Base-Base+1:Base*mp).*downch));
%     [TEraw1] = LORA.TIME_ERR(fourier1);
%     [TEraw2] = LORA.TIME_ERR(fourier2);
%     TEraw_array1(mp) = TEraw1;
%     TEraw_array2(mp) = TEraw2;
% end
% TEraw_mean1 = mean( abs(TEraw_array1) );
% TEraw_mean2 = mean( abs(TEraw_array2) );
% if(TEraw_mean1>TEraw_mean2)
%     DEraw = -1*DEraw;
% end
% for j=1:length(channel_chirp_treshold)
%     channel_chirp_treshold(j) = channel_chirp_treshold(j).*exp(1i*DEraw*j*(-1));
% end


% 1. Coarse estimation ____________________________________________________
% Определяем Hz/samp
for i = 1:num_pre
    fourier = abs(fft(channel_chirp_treshold(i*N-N+1:N*i).*downch));
    [~, ind1] = max( fourier );
%     pre_align(i) = ind1;

    if(ind1>N/2)
        pre_align(i) = (ind1-1)-N;
    else
        pre_align(i) = ind1-1;
    end
    
end


est1 = (mean(pre_align))*fps;
dphi1 = est1*2*pi*(1/BW); % сдвиг

channel_chirp_realign = zeros(1,length(channel_chirp_treshold));
for j=1:length(channel_chirp_treshold)
    channel_chirp_realign(j)=channel_chirp_treshold(j).*exp(1i*dphi1*j*(-1));
end

fourier_check = abs(fft(channel_chirp_realign(1:N).*downch));
figure(1)
stem(fourier)
title('Coarse')
% STOint
% fup_idx
% fdown_idx
% return




% 2. Дробная оценка частоты _______________________________________________
for i=1:num_pre-1
    corr_preamb = xcorr( channel_chirp_realign(N*i-N+1+N:N*i+N), channel_chirp_realign(N*i-N+1:N*i));
    [a,b]=max(corr_preamb);
    dphi_array(i) = angle(a);
end
est2 = mean(dphi_array)/(2*pi*Ts);
dphi2 = (est2*2*pi)/BW;

%     dphi2./(pi*Ts)
% est2=(dphi2(2)-dphi2(1))/(pi*Ts);
% est2=0;
% dphi2 = est2*2*pi*Ts/Base;
% dphi2=(est2*2*pi)/BW;

% return
% left_half = channel_chirp_realign(1:N/4);
% left_ref = downch(1:N/4);
% 
% right_half = channel_chirp_realign(N/4+1:N/2);
% right_ref = downch(N/4+1:N/2);
% 
% bpf3 = fft(left_half.*left_ref);
% bpf4 = fft(right_half.*right_ref);
% [max_a3] = max(bpf3);
% [max_a4] = max(bpf4);
% 
% a11=max(max_a3);
% a12=max(max_a4);
% 
% est2 = (angle(a12)-angle(a11))/(pi*Ts/2);
% dphi2 = est2*2*pi*Ts/Base;
 
% точное устранение фазового набега
for j=1:length(channel_chirp_realign)
    channel_chirp_frac_est(j) = channel_chirp_realign(j).*exp(1i*dphi2*j*(-1));
end


fourier_check = abs(fft(channel_chirp_frac_est(1:N).*downch));
figure(2)
stem(fourier_check)
title('Frac')
% return
% 3. Fine estimation ______________________________________________________
for i = 1:num_pre
    fourier2 = abs(fft(channel_chirp_frac_est(i*N-N+1:N*i).*downch));
    [~, ind1] = max( fourier2 );

    if(ind1>N/2)
        pre_align2(i) = (ind1-1)-N;
    else
        pre_align2(i) = ind1-1;
    end

end

est3 = (mean(pre_align2))*fps;
dphi3 = est3*2*pi*(1/BW); % сдвиг


            cccc = zeros(1,length(channel_chirp_frac_est));
            for j=1:length(channel_chirp_frac_est)
                cccc(j)=channel_chirp_frac_est(j).*exp(1i*dphi3*j*(-1));
            end
fourier_check = abs(fft(cccc(1:N).*downch));
figure(3)
stem(fourier_check)
title('Fine')
% corr_preamb=xcorr( channel_chirp_realign(N+1:N*2), channel_chirp_realign(1:N));
% [a,b]=max(corr_preamb);
% dphi2=angle(a);
% est2=dphi2/(2*pi*Ts);
% 
% d_phi_ocen=(est3*2*pi)/BW;
% 
% est2 = sum(est_sqr)/num_pre;
% dphi2 = est2*2*pi*Ts/Base;
% return

% est2 = 0;
% return
% est2 = (angle(a12)-angle(a11))/(pi*Ts);
% dphi2 = est2*2*pi*Ts/Base;


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
% dphi_full = dphi1+dphi2+dphi3+DEraw;
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
    
    r1 = fft(rx_preamb(1:Base).*downch);
    r11 = fft(rx_preamb(Base+1:2*Base).*downch);
    r2 = fft(rx_downch.*conj(downch));
    [~, fup_idx] = max(abs(r1));
    [~, fdown_idx] = max(abs(r2));
    if(fup_idx>(Base/2))
        fup_idx = (fup_idx-1)-Base;
    end
%     if(fdown_idx>(Base/2))
%         r2(fdown_idx)=0;
%         [~, fdown_idx] = max(abs(r2));
%     end

%     [~, fup_idx2] = max(abs(r11));
%     phi2 = angle(r1(fup_idx));
%     phi1 = angle(r2(fup_idx2));
phi2=0;
phi1=0;

end