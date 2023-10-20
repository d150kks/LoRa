clc
clear all
close all

return
tic

%% ================================= Переменные
% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
rc_size = 4;
rc = (SF-rc_size);
Base_rc = 2^rc;
rc_factor = 2^SF/Base_rc;
bits2sym = rc;
BW = 30e6;
snr = [-16:1:0];
nIter = 20;

LORA = myLoRaClass(SF,BW);
Base = LORA.Base;
downch = LORA.downch;

num_pre = 8;

num_sym = 1000;
nbits = num_sym*bits2sym; 
data = randi([0 1],1, nbits); 



% return
%% ================================= Mодуляция
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified( data, num_sym, 1);

%% ================================= BER

tic

h11 = zeros(1, length(mod_chirp));
h11(1) = 1;
h11(2) = 0.5;
h11(5) = 0.3;
H11 = fft(h11);

mod_chirp = ifft( fft(mod_chirp).*H11 );

fps = BW/Base;
freq_shift = fps*0.5;
dphi=freq_shift*2*pi*(1/BW);% сдвиг

% вводим частотный сдвиг
for j=1:length(mod_chirp)
    mod_chirp(j)=mod_chirp(j)*exp(1i*dphi*j);
end


for n = 1:length(snr)
    n 

    [numErr, NumData] = deal(0);

    for iter = 1:nIter
        
        % АБГШ 
        rxSig = awgn(mod_chirp,snr(n),'measured');

        % демодуляция
        [sv_decode, sv, fourier] = LORA.delorax_modified( rxSig, num_sym);
        dbits = int2bit(sv_decode.', rc)';
%         dbits = dbits(:)';

        % подсчет БЕР с учетом задержки
        err = sum(dbits~=data);

        % Increment the error and bit counters
        numErr = numErr + err;        
        NumData = NumData + nbits;
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

% return
% 
% save('lora_rs_ber.mat','BER')
% save('snr_crc.mat','snr')

