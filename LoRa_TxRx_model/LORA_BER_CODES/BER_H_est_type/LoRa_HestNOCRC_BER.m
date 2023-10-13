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

snr = [-16:1:0];
nIter = 100;

% массивы данных
nbits = 2100; 
num_sym = nbits/(SF);
data = randi([0 1],1, nbits); 
data = zeros(1, nbits); 

%% ================================= модуляция
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified( data, num_sym, 1);
num_pre = 4;
pre_len = num_pre*Base;
tx_preamble = repmat(LORA.chirp,1,num_pre);
tx_downch = repmat(downch,1,num_pre);
tx_chirp = [tx_downch, tx_preamble, mod_chirp]; % DOWN, UP, PAYLOAD
% tx_chirp = [downch, tx_preamble, mod_chirp]; % DOWN, UP, PAYLOAD
tx_length = length(tx_chirp);

%% ================================= BER
tic

h11 = zeros(1, Base);
h11(1) = 1;
h11(2) = 0.9;
h11(3) = 0.8;
h11(5) = 0.8;
h11(7) = 0.7;
h11(9) = 0.5;
h11(11) = 0.45;
h11(13) = 0.3;
h11(15) = 0.2;
h11(16) = 0.1;
H11 = fft(h11);
% H11 = ones(1, Base);

for i=1:tx_length/Base
    tx_sig_fft = fft(tx_chirp(i*Base-Base+1:Base*i));
    tx_chirp_h(i*Base-Base+1:Base*i) = ifft(tx_sig_fft.*H11);
end

%         N = Base;
%         for i=1:num_sym
%             fourier_data(1, i*N-N+1:i*N) = abs(fft(    tx_chirp(i*N-N+1:i*N).*1 ));
%             fourier_data(2, i*N-N+1:i*N) = abs(fft( tx_chirp_h(i*N-N+1:i*N).*1 ));
%             fourier_data(3, i*N-N+1:i*N) = abs(H11)*12;
%         end
% 
% figure(3); hold on
% stem(real(fourier_data(2,:)))
% stem(real(fourier_data(1,:)))
% plot(real(fourier_data(3,:)))
% xlim([1 Base*4])
% % return

fps = BW/Base;
freq_shift = fps*1.5;
freq_shift = fps*0.4;
dphi=freq_shift*2*pi*(1/BW);% сдвиг

% вводим частотный сдвиг
for j=1:length(tx_chirp_h)
    tx_chirp_hf(j)=tx_chirp_h(j)*exp(1i*dphi*j);
end

% BER Calculation
snr=0;
for n = 1:length(snr)

    [numErr, NumData] = deal(0);

    for iter = 1:nIter
        fprintf( 'iter num: %.d\n', n); 
        
        % АБГШ 
        delay = 2;
        rxSig = awgn([zeros(1,delay), tx_chirp_hf,zeros(1,100)], snr(n), 'measured');
        rxSig = rxSig(1:tx_length);

        [freq_data, corrected_signal, rx_preamble] = LORA.LORA_FREQ_ESTIM_v3(rxSig, num_pre);
        fprintf( 'Err:     %.2f\n', freq_shift-(freq_data{2}+freq_data{3}+freq_data{4}))
        fprintf( 'STO:     %.2f\n', freq_data{1});
        noncorrected_signal = rxSig(pre_len*2+1:end);


%         N = Base;
%         for i=1:num_sym
%             fourier_data(1, i*N-N+1:i*N) = abs(fft(    corrected_signal(i*N-N+1:i*N).*downch ));
%             fourier_data(2, i*N-N+1:i*N) = abs(fft( noncorrected_signal(i*N-N+1:i*N).*downch ));
%         end
% 
% figure(3); hold on
% stem(real(fourier_data(1,:)))
% stem(real(fourier_data(2,:)))
% return

H11hat = 0;
for i=1:num_pre
    H11hat = H11hat + fft(rx_preamble(i*Base-Base+1:i*Base))./fft(chirp);
end
H11hat = H11hat/(num_pre);

for i=1:length(corrected_signal)/Base
    symrxfft = fft(corrected_signal(i*Base-Base+1:Base*i));
%     corrected_signal(i*Base-Base+1:Base*i) = ifft( symrxfft./H11hat );
    corrected_signal(i*Base-Base+1:Base*i) = ifft( symrxfft./H11 );
end

% figure(1)
% subplot(211); hold on
% plot(real(H11))
% plot(real(H11hat))
% 
% subplot(212); hold on
% plot(imag(H11))
% plot(imag(H11hat))

% figure(2); hold on
% plot(abs(H11))
% plot(abs(H11hat))

%         N = Base;
%         for i=1:num_sym
%             fourier_data_cor(i*N-N+1:i*N) = abs(fft( corrected_signal(i*N-N+1:i*N).*1 ));
%         end
% 
% figure(4); hold on
% stem(fourier_data_cor)
% stem(real(fourier_data(1,:)))
% plot(real(fourier_data(3,:)))
% xlim([1 Base*4])

% return

        % демодуляция
        [sv_decode, sv, fourier] = LORA.delorax_modified(corrected_signal, num_sym);
        dbits = de2bi(sv_decode, SF)';
        dbits = dbits(:)';


% err_vec = (check_data~=sv_decode);
% 
%         N = Base;
%         for i=1:num_sym
%             if(err_vec(i)==1)
%                 fourier_data_cor(i*N-N+1:i*N) = abs(fft(    corrected_signal(i*N-N+1:i*N).*downch ));
%                 fourier_data_non(i*N-N+1:i*N) = abs(fft( noncorrected_signal(i*N-N+1:i*N).*downch ));
%             else
%                 fourier_data_cor(i*N-N+1:i*N) = zeros(1,Base);
%                 fourier_data_non(i*N-N+1:i*N) = zeros(1,Base);
%             end
%         end
% figure(3); hold on
% stem(real(fourier_data_non))
% stem(real(fourier_data_cor))
% % plot(real(fourier_data(3,:)))
% return
        % подсчет БЕР с учетом задержки
        err = sum(dbits~=data);

        % Increment the error and bit counters
        numErr = numErr + err;
        NumData = NumData + nbits;
        return
        clc;
    end


    % Estimate the BER for both methods
    BER(n) = numErr/NumData;
end
toc

%%
figure(1)
semilogy(snr,BER,'-*','color','k');
grid
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
hold off
title('SNR');

% save('HestNOCRC_ber.mat','BER')
% save('snr.mat','snr')

