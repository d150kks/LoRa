clc
clear all
close all


tic

%% ================================= Переменные
% Class initialization



% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
BW = 30e6;

LORA = myLoRaClass(SF, BW);
Base = LORA.Base;
downch = LORA.downch;

snr = [-16:1:0];
nIter = 100;

% массивы данных
nbits = 11200; 
num_sym = nbits/SF;
data = randi([0 1],1, nbits); 



%% ================================= модуляция
[mod_chirp, check_data] = LORA.lorax_modified( data, num_sym, 0);


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
        [sv_decode, sv, fourier] = LORA.delorax_modified(rxSig, num_sym);
        dbits = de2bi(sv,SF).';
        dbits = dbits(:).';

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
save('lora_ber.mat','BER')
save('snr.mat','snr')

