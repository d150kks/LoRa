function [sv, fourier] = delorax_modified( length_data, SF, downch, chirp)
num = length(downch);
B=2^SF;
% css = [];
for i = 1:length_data/SF

d = chirp(num*i-num+1:num*i).*downch;   % перемножаем входной и опорный ќЅ–ј“Ќый чирп


fourier = abs(fft(d));            % переводим результат в область частот
% fourier = abs(fourier(end/2:end));
[maxValue, indexMax] = max( fourier ); % находим щелчок  частоты в чирпе
% css = fourier
%%%%%%%
% вычисл€ем значение кодового слова исход€ из базы сигнала

% sv(i) = indexMax-1;
if indexMax>B
    sv(i) = B - (num - indexMax) - 1;
else
    sv(i) = indexMax - 1;
end

end

end