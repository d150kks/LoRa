function [h_1, param] = Channel(In,BW,SNR,NumOfRay,FreqShift,TimeShift)

if (nargin < 2)
    BW = 125e3;
end
if (nargin < 3)
    SNR = 0;
end
if (nargin < 4)
    NumOfRay = 10;
end
if (nargin < 5)
    FreqShift = 7500;
end
if (nargin < 6)
    TimeShift = 25;
end
%%  
param.SNR       = SNR;                                                          %  /
param.TimeShift = TimeShift;                                                    %  
param.FreqShift = FreqShift;                                                         %  
param.NumOfRay  = NumOfRay;                                                           %    
param.Ray       = [1, 0.9, 0.8, 0.8, 0.7, 0.5, 0.45, 0.3, 0.2, 0.1].*1;            %  
param.Delay     = [0,   4,   6,   7,   8,  10,   12,  15,  16,  20]+TimeShift;  %  
param.Ray(1) = 1;
%%  

dfi = param.FreqShift*2*pi/BW;
In_FreqShift = zeros(length(In),1);

for j = 1:length(In)
    In_FreqShift(j) = In(j)*exp(1i*dfi*j);
end

%%   + 

% In_TimeShift = zeros(length(In_FreqShift)+max(param.Delay),1);
% In_TimeShift(param.Delay(1):end - max(param.Delay)+param.Delay(1)-1) = In_FreqShift;
In_TimeShift = In_FreqShift;

% load("h_2.mat");
% h_time = ifft(h);

%% RELAY
NFFT    = 1024;
BW      = 2e6;
Nray    = 10;
h       = zeros(1, NFFT);
sigmaDS = 0.39;
mueDS   = -6.44;
rtau    = 2.3;
K = 10;

gI = randn(1,Nray);
gQ = randn(1,Nray);
h_2 =sqrt(1/(2*(K+1)))*(gI+1i*gQ);
h_3 = K/(K+1);
delay = -rtau*10.^(sigmaDS*randn(1,Nray)+mueDS).*log(rand(1,Nray));
delay = sort(delay-min(delay));
delay_norm = round(delay/(1/BW));
delay_uni = unique(delay_norm);

h_1 = [h_3,h_2];
% h_1 = 1;

In_TimeShift_Reley = filter(h_1,1,In_TimeShift);
% In_TimeShift = In_TimeShift_Reley(1:length(In_TimeShift));
% end




Out = awgn(In_TimeShift_Reley,SNR,'measured');
end