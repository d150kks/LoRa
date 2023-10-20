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

LORA = myLoRaClass_true(SF,BW);
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


%% ================================= CRC-RS LUT
rs_peaks_lut = (0:Base_rc-1);
rs_gray_lut = LORA.grayCode(rs_peaks_lut+1);
rs_int_lut = int2bit(rs_gray_lut.', rc).';
rs_crc_lut = LORA.codeCRC(rs_int_lut, Base_rc);
rs_int_crc_lut = bit2int(rs_crc_lut.', SF).';
rs_lut = rs_int_crc_lut(rs_gray_lut+1);


%% ================================= CRC
data_crc = LORA.codeCRC(data, num_sym);


%% ================================= Mодуляция
% формирование пустых массивов
mod_chirp = zeros(1, num_sym*Base);           % чирпы

% Модуляция
for i = 1:num_sym
    code_word = bit2int(data_crc(SF*i-SF+1:SF*i).', SF);      % значение кодового слова
    code_word_rs = (find(rs_lut==code_word)-1)*rc_factor;
    check_data(i) = code_word_rs;
    

    cs = single((code_word_rs/Base)*Ts/ts);             % место сдвига
    chirp1 = chirp(cs+1:end); % первый чирп
    chirp2 = chirp(1:cs); % второй чирп
    mod_chirp(i*Base-Base+1:Base*i) = [chirp1, chirp2];
end 


%% ================================= Демодуляция

for i=1:num_sym
    fourier(i,:) = abs(fft(mod_chirp(Base*i-Base+1:Base*i).*downch));   
    fourier_rs = LORA.reduced_set_fourier(fourier(i,:));
    [~, indexMax] = max( fourier_rs ); % находим щелчок  частоты в чирпе
    sv(i) = rs_lut(indexMax);

%     [sv, sv_cor, peakMakcor, dbits] = HARD_CRC_DEMOD( fourier_rs, LORA );
%     sv1(i) = sv;
% return
end
% dbits = int2bit(sv.', rc)';
% err = sum(dbits~=data)

check_data
sv
sv_bits = int2bit(sv.', SF).';

% [fourier_rc] = reduced_set_fourier(fourier, Base_rc, rc_factor, LORA);
% [fourier_rc] = LORA.reduced_set_fourier(fourier);
%
figure(1); hold on
stem(fourier.')
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


i=1;
fourier = abs(fft(rxSig(Base*i-Base+1:Base*i).*downch));   

% [fourier_rc] = reduced_set_fourier(fourier, Base_rc, rc_factor, LORA);
[fourier_rc] = LORA.reduced_set_fourier(fourier);

figure(1); hold on
stem(fourier)
stem( fourier_rc)
% stem(fourier_cconv(1,:))
% stem(fourier.')
% stem(fourier_cconv.');

return


% демодуляция

        % Param
        Gcode = LORA.grayCode;
        hard_bits = zeros(1,SF*num_sym);
        
        sv = zeros(1,num_sym);
        sv_cor = zeros(1,num_sym);
        soft_bits = [];
        



%%
figure(1)
semilogy(snr,BER,'-*','color','k');
grid
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
hold off
title('SNR');



function [fourier_rc] = reduced_set_fourier(fourier, Base_rc, rc_factor, LORA)
    
    % Reduced set vectors and aos
    rc_peaks = (0:Base_rc-1)*rc_factor+1;
    rc_aos = -rc_factor/2:rc_factor/2;
    fourier_rc = zeros(1,Base_rc);

    for rc_peak_idx=1:Base_rc
        rc_win = LORA.CYC_SHIFT(rc_peaks(rc_peak_idx)+rc_aos);
        fourier_rc(rc_peak_idx) = sum(fourier(rc_win))/sqrt(rc_factor);
    end

end

