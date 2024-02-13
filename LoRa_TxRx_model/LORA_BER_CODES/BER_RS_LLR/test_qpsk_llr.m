clc
clear all
close all

data = randi([0, 1],1 , 40);
data_mod = qammod(data.',4, InputType="bit").';

[data_rx, nvar] = awgn(data_mod, 20, 'measured');

data_demod = -qamdemod(data_rx.', 4, OutputType='approxllr', NoiseVariance=nvar).';

figure(1); hold on
plot( normalize(data) )
plot( normalize(data_demod) )

