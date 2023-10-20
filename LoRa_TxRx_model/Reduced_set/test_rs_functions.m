clc
clear all
close all


tic

%% ================================= Переменные
% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
rc_size = 4;
rc = (SF-rc_size);
Base_rc = 2^rc;
rc_factor = 2^SF/Base_rc;
b2s = rc;

BW = 125e3;
snr = 5;

LORA = myLoRaClass(SF,BW);
Base = LORA.Base;
downch = LORA.downch;
chirp = LORA.chirp;
Ts = LORA.Ts;
ts = LORA.ts;

num_pre = 8;
num_sym = Base_rc;
% nbits = b2s*num_sym; 
% data = randi([0 1],1, nbits); 
datade = 0:Base_rc-1;
data = int2bit(datade.', rc).';



%% ================================= Mодуляция
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified(data, num_sym, 1);
[sv_decode, sv, fourier] = LORA.delorax_modified( mod_chirp, num_sym);

return
%% ================================= Канал

h11 = zeros(1, length(mod_chirp));
h11(1) = 1;
h11(2) = 0.5;
h11(5) = 0.3;
H11 = fft(h11);

mod_chirp = ifft( fft(mod_chirp).*H11 );

fps = BW/Base;
freq_shift = fps*1.5;
dphi=freq_shift*2*pi*(1/BW);% сдвиг

% вводим частотный сдвиг
for j=1:length(mod_chirp)
    mod_chirp(j)=mod_chirp(j)*exp(1i*dphi*j);
end

rxSig = awgn(mod_chirp, snr, 'measured');


%% ================================= Демодуляция
% 1. type demod
% [sv_decode, sv, fourier] = LORA.delorax_modified(mod_chirp, num_sym);
% sum(abs(sv_decode-check_no_gray))
% sum(abs(sv-check_data/4))

% 2. type demod
[soft_bits, hard_bits, sv, sv_cor, fourier] = LORA.DELORAX_CRC(mod_chirp, num_sym, 1, 2);
sum(abs(data-hard_bits))


% i=1;
% [sv_decode, sv, fourier] = LORA.delorax_modified(mod_chirp(Base*i-Base+1:Base*i), 1);

% figure(1); hold on
% plot(data)
% plot(soft_bits)

return

return


