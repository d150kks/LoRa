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
% data = randi([0 1],1, Base_rc*rc); 

%% ================================= CRC-RS LUT
rs_peaks_lut = (0:Base_rc-1);
rs_gray_lut = LORA.grayCode(rs_peaks_lut+1);



%% ================================= CRC
% num_sym = 1;
% data = [1 0 1];
% data = [0 0 1 0 0 0 0];
% data_crc = LORA.codeCRC(data, num_sym);

% data_crc_bad = [1 0 1 1 1 1 0]
% sum(LORA.CRC4(data_crc_bad.').')
% return

%% ================================= Mодуляция
% формирование пустых массивов
mod_chirp = zeros(1, num_sym*Base);           % чирпы

% Модуляция
for i = 1:num_sym
    code_word = bit2int(data(rc*i-rc+1:rc*i).', rc);      % значение кодового слова

    code_word_rs = (rs_peaks_lut(code_word==rs_gray_lut))*rc_factor;
    check_data(i) = code_word_rs;
    code_word_data(i) = code_word;

    cs = single((code_word_rs/Base)*Ts/ts);             % место сдвига
    chirp1 = chirp(cs+1:end); % первый чирп
    chirp2 = chirp(1:cs); % второй чирп
    mod_chirp(i*Base-Base+1:Base*i) = [chirp1, chirp2];

end 
% check_data
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified_crcrs(data, num_sym);
% mod_chirp = circshift(mod_chirp, 8);
check_data
delay1 = 32;
delay2 = 68;
% mod_chirp = mod_chirp + [zeros(1,delay1), mod_chirp(1:end-delay1)] + [zeros(1,delay2), mod_chirp(1:end-delay2)];
mod_chirp = awgn(mod_chirp, -0, 'measured', 17);

%% ================================= Демодуляция
for i=1:num_sym
    d = [-1, 1, -1];
%     d = gausswin(9).';
%     d(1:2:end)=0;
%     d = repmat(d,1,64);
    fourier(i,:) = abs(fft(mod_chirp(Base*i-Base+1:Base*i).*downch)); 
    fourier_rs = abs(filtfilt( d/8, 1, fourier(i,:).*1 ));
    fourier_rs = 128*fourier_rs./max(fourier_rs);
    fourier_rs2 = LORA.reduced_set_fourier( fourier_rs );

%     fourier_rs = LORA.reduced_set_fourier(fourier(i,:));
    [~, indexMax] = max( fourier_rs2 ); % находим щелчок  частоты в чирпе
    sv(i) = LORA.grayCode(indexMax);
    sv_rs(i) = sv(i)*rc_factor;

end
% [sv_rs, sv, fourier] = LORA.delorax_modified( mod_chirp, num_sym);
% [hard_bits, sv_rs, sv, fourier, fourier_rs] = LORA.delorax_crcrs( mod_chirp, num_sym);
%

% check_data;
% sv_rs
% [data_decode] = LORA.decodeCRC(sv_rs_crc, num_sym, 0, 0);
hard_bits = int2bit(sv.', rc).';
err = sum(hard_bits~=data)



figure(1); hold on
% plot( real(wt) )
stem( fourier_rs, 'b')
stem( fourier(i,:).', 'r')
return



function [y] = my_fft(x)
    N = length(x);
    y = zeros(1,N);
    for k=1:N
        for j=1:N
            y(k) = y(k)+x(j).*exp(-1i*2*pi*k*j/N);
        end
    end
end

function [y] = my_chirplet(x)
    N = length(x);
    y = zeros(1,N);
    m=1/(N^1);
    for k=1:N
        for j=1:N
            y(k) = y(k)+x(j).*exp(-1i*2*pi*k*m*j.^2/2 );
        end
    end
end