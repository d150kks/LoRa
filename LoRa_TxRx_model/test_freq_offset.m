clc
clear all
close all


tic

%% ================================= Переменные
% Class initialization



% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
BW = 125e3;

LORA = myLoRaClass(SF, BW);
Base = LORA.Base;
downch = LORA.downch;
chirp = LORA.chirp;
Ts = LORA.Ts;

nIter = 10;

% массивы данных
nbits = 1120; 
num_sym = nbits/SF;
data = randi([0 1],1, nbits); 
data = zeros(1,nbits);


%% ================================= Rate matching
[dataRM, numcodebitsRM, num_sym, zeros2end, flagRM] = LORA.RM(data);


%% ================================= CRC coding
[data_code] = LORA.codeCRC(dataRM, num_sym);

%% ================================= модуляция
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified( data_code, num_sym, 1);
num_pre = 4;
pre_len = Base*num_pre;
tx_preamble = repmat(chirp,1,num_pre);
tx_downch = repmat(downch,1,num_pre);
tx_chirp = [tx_downch, tx_preamble, mod_chirp]; % DOWN, UP, PAYLOAD
tx_length = length(tx_chirp);
% sig_length = tx_length-Base;

%% ================================= Channel
tic

% ~~~~~~~~ Freq multipath ~~~~~~~~ 
h11 = zeros(1, Base);
% h11 = zeros(1, tx_length);
h11(1) = 1;
h11(2) = 0.5;
h11(5) = 0.3;
H11 = fft(h11);
% H11 = ones(1, Base);

% h11 = zeros(1, Base);
% h11(1) = 1;
% h11(2) = 0.9;
% h11(3) = 0.8;
% h11(5) = 0.8;
% h11(7) = 0.7;
% h11(9) = 0.5;
% h11(11) = 0.45;
% h11(13) = 0.3;
% h11(15) = 0.2;
% h11(16) = 0.1;
% H11 = fft(h11);
% H11 = ones(1, Base);

for i=1:tx_length/Base
    tx_sig_fft = fft(tx_chirp(i*Base-Base+1:Base*i));
    tx_chirp_h(i*Base-Base+1:Base*i) = ifft(tx_sig_fft.*H11);
end

% tx_chirp = ifft( fft(tx_chirp).*H11 );
% tx_chirp = tx_chirp;

% ~~~~~~~~ Freq Shift ~~~~~~~~ 
fps = BW/Base;
% freq_shift = fps*1.5;
% freq_shift = fps*1.5;
% freq_shift = fps*3.5;
freq_shift = fps*1.3;
% freq_shift = fps*0.5;
freq_shift = 0;
dphi=freq_shift*2*pi*(1/BW);% сдвиг

% вводим частотный сдвиг
for j=1:length(tx_chirp_h)
    tx_chirp_hf(j)=tx_chirp_h(j)*exp(1i*dphi*j);
end

% ~~~~~~~~ Awgn ~~~~~~~~ 
% for pp=1:100;
delay = 4*1;
snr = 30;
% pq = 100;
% tx_chirp_hf = resample(tx_chirp_hf,pq,1);
rx_sig_hfn = awgn( [zeros(1,delay), tx_chirp_hf, zeros(1,Base)], snr, 'measured');
% rx_sig_hfn = rx_sig_hfn(65:pq:end);

%% ================================= Receiver
rx_downch = rx_sig_hfn(1:pre_len);
rx_pre = rx_sig_hfn(pre_len+1:pre_len*2);
rx_payload = rx_sig_hfn(pre_len*2+1:end);
%
%
%
%
% ~~~~~~~~ 0. Coarse CFO and STO estimation ~~~~~~~~ 
N = Base;
OS = 10;
r1 = 0;
r2 = 0;
for i = 1:num_pre
    r1 = r1 + abs( fft([rx_pre(i*N-N+1:N*i).*downch, zeros(1,N*(OS-1))]) );
    r2 = r2 + abs( fft([rx_downch(i*N-N+1:N*i).*chirp, zeros(1,Base*(OS-1))]) );
end

% r2 = fft( [rx_downch.*conj(downch), zeros(1,Base*(OS-1))] );
[~, fup_idx] = max(abs(r1));
[~, fdown_idx] = max(abs(r2));

if(fup_idx>(Base*OS/2))
    fup_idx = (fup_idx-1)-Base*OS;
end
if(fdown_idx>(Base*OS/2))
    fdown_idx = (fdown_idx-1)-Base*OS;
end
fup_idx = fup_idx/OS; 
fdown_idx = fdown_idx/OS;

STO = (fup_idx-fdown_idx)/2;
STOint = round(STO);
STOfraq = STO-STOint;
CFO = fps*(fup_idx+fdown_idx)/2;
CFOdphi = CFO*2*pi*(1/BW);

% figure(1); hold on
% plot(abs(r1))
% plot(abs(r2))
% return

% ~~~~~~~~ 0. Coarse CFO and STO compensation ~~~~~~~~ 
% rx_sig_hfn2 = rx_sig_hfn;
rx_sig_hfn_fft = fft(rx_sig_hfn);
STOintphi = STO*2*pi/length(rx_sig_hfn_fft);
% STOintphi = STOint*2*pi/length(rx_sig_hfn_fft);
for j=1:length(rx_sig_hfn_fft)
    rx_sig_hfn_fft(j) = rx_sig_hfn_fft(j).*exp(1i*j*-STOintphi);
end
rx_sig_hfn2 = ifft(rx_sig_hfn_fft);
% rx_sig_hfn2 = LORA.DPHI_COMP(rx_sig_hfn2, CFOdphi);


% figure(1)
% plot(real(tx_chirp))
% hold on
% plot(real(rx_sig_hfn2))
% xlim([1, Base*1])
% return

rx_sig_hfn2 = LORA.DPHI_COMP(rx_sig_hfn2, CFOdphi);
rx_downch2 = rx_sig_hfn2(1:pre_len);
rx_pre2 = rx_sig_hfn2(pre_len+1:pre_len*2);
rx_payload2 = rx_sig_hfn2(pre_len*2+1:end);


% ~~~~~~~~ 2. Fraq estimation ~~~~~~~~ 
% По двум последовательным одинаковым чирпам как в (1.) вычисляем CFO_fraq
left_ref  = downch(1:N/2);
right_ref = downch(N/2+1:N);
for i = 1:num_pre
    chirp2est = rx_pre2(i*N-N+1:i*N);
    left_half = chirp2est(1:N/2);
    right_half = chirp2est(N/2+1:N);
    
    bpf3 = fft(left_half.*left_ref);
    bpf4 = fft(right_half.*right_ref);

    [~, max_a3] = max(abs(bpf3));
    [~, max_a4] = max(abs(bpf4));

    [a11] = bpf3(max_a3);
    [a12] = bpf4(max_a4);

    est2_reg(i) = (angle(a12)-angle(a11));
end
est2_reg = 0;
est2 = mean(est2_reg)/(pi*Ts);
dphi2 = est2*2*pi*Ts/Base;
rx_sig_hfn2 = LORA.DPHI_COMP(rx_sig_hfn2, dphi2);


% ~~~~~~~~ 3. fine estimation ~~~~~~~~ 
% Устранение фазового сдвига
aos = -2:2;
argumon = zeros(1, num_pre-1);
for i = 1:num_pre-1
    Y1 = fft(rx_pre2(i*N+1:N*i+N).*downch);
    Y0 = fft(rx_pre2(i*N-N+1:N*i).*downch);
    [ampmax1, indmax1] = max(Y1);
    [ampmax0, indmax0] = max(Y0);
    indwin1 = LORA.CYC_SHIFT(indmax1+aos);
    indwin0 = LORA.CYC_SHIFT(indmax0+aos);

    argumon(i) = angle( sum( Y1(indwin1) .* conj(Y0(indwin0)) ) );
end
% argumon
[arg] = mean(argumon);
est3 = arg/(2*pi*Ts);
dphi3 = est3*2*pi/BW; % сдвиг
% ~~~~~~~~ 3. fine Compensation ~~~~~~~~ 
rx_sig_hfn3 = LORA.DPHI_COMP(rx_sig_hfn2, dphi3);
rx_downch3 = rx_sig_hfn3(1:Base);
rx_pre3 = rx_sig_hfn3(Base+1:Base*9);
rx_payload3 = rx_sig_hfn3(Base*9+1:end);


H11hat = 0;
for i=1:num_pre
%     H11hat = H11hat + fft(rx_pre(i*Base-Base+1:i*Base))./fft(chirp);
    H11hat = H11hat + fft(rx_pre3(i*Base-Base+1:i*Base))./fft(chirp);
end
H11hat = H11hat/(num_pre);

figure(1)
subplot(211); hold on
plot(real(H11))
plot(real(H11hat))

subplot(212); hold on
plot(imag(H11))
plot(imag(H11hat))

figure(2); hold on
plot(abs(H11))
plot(abs(H11hat))


fprintf( 'Offset:  %.2f\n', freq_shift)
fprintf(['STO:     %.2f\n' ...
         'STOint:  %.2f\n' ...
         'STOfraq: %.2f\n'], STO, STOint, STOfraq)
fprintf( 'Est CFO: %.2f\n', CFO)
fprintf( 'Est Frq: %.2f\n', est2)
fprintf( 'Est Fin: %.2f\n', est3)
fprintf( 'Err:     %.2f\n\n', freq_shift-(CFO+est2+est3))

% errvec(pp) = freq_shift-(CFO+est3);
% if(abs(errvec(pp))>50)
%     break
% end
% end
% 
% plot(errvec)
% return

for i=1:length(rx_payload3)/Base
    symrxfft = fft(rx_payload3(i*Base-Base+1:Base*i));
%     corrected_signal(i*Base-Base+1:Base*i) = ifft( symrxfft./H11hat );
    corrected_signal(i*Base-Base+1:Base*i) = ifft( symrxfft./H11 );
end


% % ~~~~~~~~ 2. Fraq estimation ~~~~~~~~ 
%             % По двум последовательным одинаковым чирпам как в (1.) вычисляем CFO_fraq
%             left_ref  = downch(1:N/2);
%             right_ref = downch(N/2+1:N);
%             est2_reg = zeros(1,num_pre);
%             for i = 1:num_pre
%                 chirp2est = rx_pre(i*N-N+1:i*N);
%                 left_half = chirp2est(1:N/2);
%                 right_half = chirp2est(N/2+1:N);
% 
%                 bpf3 = fft( left_half.*left_ref);
%                 bpf4 = fft( right_half.*right_ref);
%                 
%                 [max_a3] = max(bpf3);
%                 [max_a4] = max(bpf4);
%                 
%                 a11=max(max_a3);
%                 a12=max(max_a4);
%                 est2_reg(i) = (angle(a12)-angle(a11));
%             end
% %             est2_reg
%             est2 = mean(est2_reg)/(pi*Ts);
%             dphi2 = est2*2*pi*Ts/Base;
% % rx_pre2 = LORA.DPHI_COMP(rx_pre, dphi2);
% rx_pre2 = rx_pre;
% % ~~~~~~~~ 3. fine estimation ~~~~~~~~ 
% % Устранение фазового сдвига
% N = Base;
% aos = -2:2;
% argumon = zeros(1, num_pre-1);
% for i = 1:num_pre-1
%     Y1 = fft(rx_pre2(i*N+1:N*i+N).*downch);
%     Y0 = fft(rx_pre2(i*N-N+1:N*i).*downch);
%     [ampmax1, indmax1] = max(Y1);
%     [ampmax0, indmax0] = max(Y0);
%     indwin1 = LORA.CYC_SHIFT(indmax1+aos);
%     indwin0 = LORA.CYC_SHIFT(indmax0+aos);
% 
%     argumon(i) = angle( sum( Y1(indwin1) .* conj(Y0(indwin0)) ) );
% end
% pos = abs(sum(argumon>=0));
% neg = abs(sum(argumon<0));
% if( pos>=neg )
%   argumon_update = argumon(argumon>=0);
% else
%   argumon_update = argumon(argumon<0);
% end
% plot(argumon_update)
% [arg] = mean(argumon_update);
% est3 = arg/(2*pi*Ts);
% dphi3 = est3*2*pi/BW; % сдвиг


% демодуляция
[soft_bits, hard_bits, sv, sv_cor, fourier] = LORA.DELORAX_CRC(corrected_signal, num_sym, tx_preamble, rx_pre3);

% CRC decoding
[data_decode] = LORA.decodeCRC(hard_bits, num_sym, zeros2end, flagRM);

% подсчет БЕР с учетом задержки
err = sum(data_decode~=data);
err
