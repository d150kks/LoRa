clc
clear all
close all


% N = 10000;
fs = 10e6;
dec_bit = 50;


f0 = 5e6;
f2 = f0+5e6;
f1 = f0-5e6;
% f3 = f1+5;
BW = f2-f1;
sf = 12;
N = fs*2^sf/BW;
T = N/fs;
t = (0:N-1)/fs;

% sf=log2(BW*N/fs);


mu = 2*pi*BW/T;

chirp = sin(2*pi*f1*t+mu*t.^2/2);
% chirp3 = sin(2*pi*f3*t+mu*t.^2/2);

figure(1)
plot(abs(fft(chirp)))
return
time_val = 1/fs:1/fs:N/fs;
figure
plot(time_val,((chirp)));
set(gca,'FontName','Times New Roman','fontsize',12);
grid on
ylabel('Амплитуда');
xlabel('Время, с');
%
sp = fft(chirp);
%
shag_freq = fs/N;
freq_val = shag_freq:shag_freq:fs;
figure
plot(freq_val, abs(sp));
set(gca,'FontName','Times New Roman','fontsize',12);
grid on
ylabel('Амплитуда');
xlabel('Частота, Гц');


shag_freq1 = f2/N;
freq_val1 = shag_freq1:shag_freq1:f2;

figure
plot(time_val,freq_val1);

set(gca,'FontName','Times New Roman','fontsize',12);
grid on

ylabel('Частота, Гц');
xlabel('Время, с');



% fshift = BW*95/(2^sf)

freq_step_1 =  (f2-f1)/N;



% формируем кодовые слова



n_bit = sf;%(2^sf)*sf;
% bits = [1 1 1 0 1 0 1 0 0 1];%randint (1,n_bit);
%bi2de(bits)

% формируем массив временных отрезков кратных sf

n_start = (N/2^sf)*dec_bit;

freq_znach  = n_start*freq_step_1;

A = 1:20;

chirp_shift = circshift(chirp,[0,-n_start+N/2^sf]);




% freq_val1_shift = circshift(freq_val1,[0,-n_start+N/2^sf]);
%
% figure
% plot(time_val,freq_val1_shift);
%
% set(gca,'FontName','Times New Roman','fontsize',12);
% grid on
%
% ylabel('Частота, Гц');
% xlabel('Время, с');
%
%
%
%
% figure
% plot(time_val,((chirp_shift)));
% set(gca,'FontName','Times New Roman','fontsize',12);
% grid on
%
% ylabel('Амплитуда');
% xlabel('Время, с');


%% канал


chirp_shift = circshift(chirp,[0,-n_start+N/2^sf]);
t_sdvig = 1e-99;
 chirp_shift0 = time_sdvig(chirp_shift, t_sdvig, fs);

 t_sdvig = 5e-8;
 chirp_shift1 = time_sdvig(chirp_shift, t_sdvig, fs)*0.8;
 t_sdvig = 1e-6;
 chirp_shift2 = time_sdvig(chirp_shift, t_sdvig, fs)*0.2;


% chirp_shift0 = [chirp_shift, zeros(1,6)];
% chirp_shift1 = [zeros(1,3),chirp_shift,zeros(1,3)]*0.7;
% chirp_shift2 = [zeros(1,6),chirp_shift]*0.3;

chirp_shift_channel = chirp_shift0+chirp_shift1+chirp_shift2;
chirp_shift_channel = chirp_shift_channel(1:length(chirp_shift));
chirp_shift_channel = awgn(chirp_shift_channel,-17,'measured');

chirp_noise = awgn(chirp,-10,'measured');
korr = xcorr(chirp,chirp_noise);

figure
plot(abs(korr));



 chirp_inv  = chirp(end:-1:1);
  chirp_pr  =(chirp_shift_channel.*chirp_inv);


%   dw =
%   chirp_pr_1 = cos(


 chirp_pr_sp  =fft(chirp_pr);



 figure
plot(time_val,((chirp_pr)));
set(gca,'FontName','Times New Roman','fontsize',12);
grid on

ylabel('Амплитуда');
xlabel('Время, с');
%
%
  figure
 plot((abs(chirp_pr_sp)));
 grid on
set(gca,'FontName','Times New Roman','fontsize',12);
grid on

ylabel('Амплитуда');
xlabel('Частота, Гц');
%


[AA,BB] = max(abs(chirp_pr_sp(1:N/2)));
if BB > 2^sf
BB = BB-2^sf;
end
BB = BB+1;

%
 demod_bi = de2bi(BB);
 demod_bit = zeros(1,sf);
 demod_bit(1:length(demod_bi)) = demod_bi;

 spec_eff = N/sf

%
%
%



