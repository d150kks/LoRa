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
rs_int_lut = int2bit(rs_gray_lut.', rc).';
rs_crc_lut = LORA.codeCRC(rs_int_lut, Base_rc);
rs_int_crc_lut = bit2int(rs_crc_lut.', SF).'
rs_lut = rs_peaks_lut*rc_factor;


%% ================================= CRC
% num_sym = 1;
% data = [1 0 1];
% data = [0 0 1 0 0 0 0];
data_crc = LORA.codeCRC(data, num_sym);

% data_crc_bad = [1 0 1 1 1 1 0]
% sum(LORA.CRC4(data_crc_bad.').')
% return

%% ================================= Mодуляция
% формирование пустых массивов
mod_chirp = zeros(1, num_sym*Base);           % чирпы

% Модуляция
for i = 1:num_sym
    code_word = bit2int(data_crc(SF*i-SF+1:SF*i).', SF);      % значение кодового слова
    code_word_rs = rs_lut(rs_int_crc_lut==code_word);
    check_data(i) = code_word_rs;

    cs = single((code_word_rs/Base)*Ts/ts);             % место сдвига
    chirp1 = chirp(cs+1:end); % первый чирп
    chirp2 = chirp(1:cs); % второй чирп
    mod_chirp(i*Base-Base+1:Base*i) = [chirp1, chirp2];

end 
% check_data
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified_crcrs(data_crc, num_sym);
mod_chirp = circshift(mod_chirp, 8);
% check_data
% return

%% ================================= Демодуляция
for i=1:num_sym
    fourier(i,:) = abs(fft(mod_chirp(Base*i-Base+1:Base*i).*downch));   
    fourier_rs = LORA.reduced_set_fourier(fourier(i,:));
    [~, indexMax] = max( fourier(i,:) ); % находим щелчок  частоты в чирпе
%     sv_rs(i) = ( LORA.grayCode(indexMax)+1 );

    [sv, sv_cor, peakMakcor, dbits] = HARD_CRC_DEMOD( fourier(i,:), LORA, rs_int_crc_lut);
    sv_cor
    sv_rs(i) = rs_lut(sv_cor==rs_int_crc_lut);
% return
end
% [sv_rs, ~, ~] = LORA.delorax_crcrs( mod_chirp, num_sym );


% check_data;
% sv_rs
sv_rs_crc = int2bit(sv_rs', SF)';
[data_decode] = LORA.decodeCRC(sv_rs_crc, num_sym, 0, 0);
% sv_bits = int2bit(sv_rs.', rc).';
% err = sum(data_decode~=data)

% [fourier_rc] = reduced_set_fourier(fourier, Base_rc, rc_factor, LORA);
% [fourier_rc] = LORA.reduced_set_fourier(fourier);
%
figure(1); hold on
% stem(fourier.')
stem(fourier_rs)



function [sv, sv_cor, peakMakcor, dbits] = HARD_CRC_DEMOD( fourier_rs, LORA, rs_int_crc_lut)
aos = 1;
aos_win = -aos:aos;
Base = 128;
grayCode = LORA.grayCode;
SF = 7;

    % ~~~~~~~~ Initial conditions ~~~~~~~~
    check_crc = 1;

    % ~~~~~~~~ CRC Demodulation ~~~~~~~~
    while check_crc~=0
    
        % Circular shift in peak search
        [~, indexMax] = max( fourier_rs ); % находим щелчок  частоты в чирпе
        sv = indexMax;

        
        peak_win = indexMax+aos_win;

        for pk=1:length(peak_win)
            if(peak_win(pk)<=0)
                peak_win(pk)=peak_win(pk)+Base;
            end
            if(peak_win(pk)>Base)
                peak_win(pk)=peak_win(pk)-Base;
            end
        end


        % Set type of the peak sort
        peaks_amp = fourier_rs(peak_win);
        peak_sort = zeros(1,2*aos+1);
        peak_sort(1) = indexMax;

        for pk=1:aos
            if( peaks_amp(aos+1+pk)>peaks_amp(aos+1-pk) )
                peak_sort(2*pk)       = peak_win(aos+1+pk);
                peak_sort(2*(pk+1)-1) = peak_win(aos+1-pk);
            else
                peak_sort(2*pk)       = peak_win(aos+1-pk);
                peak_sort(2*(pk+1)-1) = peak_win(aos+1+pk);
            end
        end

        sort_amp = fourier_rs(peak_sort);
        peak_sort(sort_amp==0)=[];
        sort_amp(sort_amp==0)=[];

% dbits
% check_crc


        % Search peak according to the crc
        for n=1:length(peak_sort)
%             sv_cor = find(grayCode==peak_sort(n)-1)-1
            sv_cor = rs_int_crc_lut(LORA.grayCode(peak_sort(n))+1);

sv_cor
            dbits = int2bit(sv_cor', SF)';

            check_crc = sum(LORA.CRC4(dbits.').');

% sv_cor
% check_crc
% return
            if(check_crc==0)
                peakMakcor = sort_amp(n);
                break
            else
                fourier_rs(peak_sort(n))=0;
            end
        end
    end
end