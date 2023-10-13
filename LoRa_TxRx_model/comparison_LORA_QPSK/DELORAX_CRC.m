function [demod_bits, sv_cor, sv, fourier] = DELORAX_CRC( length_data, SF, downch, chirp, aos)
    num = length(downch);
    B=2^SF;
    demod_bits = [];
    aos_win = -aos:aos;

    for i = 1:length_data/SF

        % Fourier
        d = chirp(num*i-num+1:num*i).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
        fourier = abs(fft(d));            % переводим результат в область частот
        % fourier = abs(fourier(end/2:end));
        [~, indexMax] = max( fourier ); % находим щелчок  частоты в чирпе

        % вычисляем значение кодового слова исходя из базы сигнала
        
        % CRC
        sv(i) = indexMax-1;
        peak_win = indexMax+aos_win;
        peak_win(peak_win<=0)=[];
        peak_win(peak_win>B)=[];
        peaks_amp = fourier(peak_win);
        [~, sort_idx] = sort(peaks_amp,'descend');
        peak_sort = peak_win(sort_idx);
        for n=1:length(peak_sort)
            sv_cor(i) = peak_sort(n)-1;
            dbits = de2bi(sv_cor(i),SF)';
            dbits = dbits(:)';
            check_crc = sum(CRC4( dbits.').');
            if(check_crc==0)
                break
            end
        end
        demod_bits = [demod_bits, dbits];
    end

end