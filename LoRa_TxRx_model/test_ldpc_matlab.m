clc
clear all
close all



%% ================================= Переменные
NumInformationBits = 486;
infoBits = rand(NumInformationBits,1) < 0.5;
P = [16 17 22 24  9  3 14 -1  4  2  7 -1 26 -1  2 -1 21 -1  1  0 -1 -1 -1 -1
     25 12 12  3  3 26  6 21 -1 15 22 -1 15 -1  4 -1 -1 16 -1  0  0 -1 -1 -1
     25 18 26 16 22 23  9 -1  0 -1  4 -1  4 -1  8 23 11 -1 -1 -1  0  0 -1 -1
      9  7  0  1 17 -1 -1  7  3 -1  3 23 -1 16 -1 -1 21 -1  0 -1 -1  0  0 -1
     24  5 26  7  1 -1 -1 15 24 15 -1  8 -1 13 -1 13 -1 11 -1 -1 -1 -1  0  0
      2  2 19 14 24  1 15 19 -1 21 -1  2 -1 24 -1  3 -1  2  1 -1 -1 -1 -1  0
    ];

blockSize = 27;
pcmatrix = ldpcQuasiCyclicMatrix(blockSize,P);



%% ================================= LDPC coding
cfgLDPCEnc = ldpcEncoderConfig(pcmatrix);
codeword = ldpcEncode(infoBits, cfgLDPCEnc);

% codeword(10:12) = ~codeword(10:12);
% codeword(10) = ~codeword(10);
% codeword(131) = ~codeword(131);
% codeword(256) = ~codeword(256);


%% ================================= LDPC decoding
cfgLDPCDec = ldpcDecoderConfig(pcmatrix);


M = 4;
maxnumiter = 10;
snr = 5;
numframes = 10;

data = randi([0 1],NumInformationBits,1,'int8');
% Transmit and receive with LDPC coding
data_bad = data;
encodedData = ldpcEncode(data_bad,cfgLDPCEnc);
encodedDatabad = encodedData;
encodedDatabad(10) = ~encodedDatabad(10);
encodedDatabad(131) = ~encodedDatabad(131);
encodedDatabad(256) = ~encodedDatabad(256);

modSignal = pskmod(data,M,InputType='bit');
% [rxsig, noisevar] = awgn(modSignal,snr);
rxsig = modSignal;
noisevar = 0.1;
demodSignal = pskdemod(rxsig,M, ...
            OutputType='approxllr', ...
            NoiseVariance=noisevar);
demodSignal2 = pskdemod(rxsig,M, ...
            OutputType='approxllr', ...
            NoiseVariance=10);
% rxbits = ldpcDecode(demodSignal,cfgLDPCDec,maxnumiter);
% errStats = sum(data~=rxbits);
% 
% fprintf(['SNR = %2d\n' ...
%         'Number of errors = %d\n'], ...
%         snr,errStats)


% figure(1)
% plot((demodSignal))
% hold on
% plot((demodSignal2))
% return
%%%%%%
% bps = 2;
S = [0,0;...
     0,1;...
     1,0;...
     1,1];

S0 = [1+1i*0, 0+1i, 0-1i, -1+1i*0];

for i=1:2
    idx_s0 = find(S(:,i)==0);
    idx_s1 = find(S(:,i)==1);
    sx0(i,:) = real(S0(idx_s0));
    sy0(i,:) = imag(S0(idx_s0));
    sx1(i,:) = real(S0(idx_s1));
    sy1(i,:) = imag(S0(idx_s1));
end

% sx01 = [real(S0(1)), real(S0(2))];
% sy01 = [imag(S0(1)), imag(S0(2))];
% sx11 = [real(S0(3)), real(S0(4))];
% sy11 = [imag(S0(3)), imag(S0(4))];
% 
% sx02 = [real(S0(1)), real(S0(3))];
% sy02 = [imag(S0(1)), imag(S0(3))];
% sx12 = [real(S0(2)), real(S0(4))];
% sy12 = [imag(S0(2)), imag(S0(4))];

L = [];
for b=1:NumInformationBits/2
    x = real(rxsig(b));
    y = imag(rxsig(b));

    for i=1:2
%     LLR1 = log( exp(-( (x-sx01).^2+(y-sy01).^2 ))/exp(-( (x-sx11).^2+(y-sy11).^2 )) );
%     LLR2 = log( exp(-( (x-sx02).^2+(y-sy02).^2 ))/exp(-( (x-sx12).^2+(y-sy12).^2 )) );
%     LLR1 =  -(min( (x-sx01).^2+(y-sy01).^2 ) - min( (x-sx11).^2+(y-sy11).^2 ));
%     LLR2 =  -(min( (x-sx02).^2+(y-sy02).^2 ) - min( (x-sx12).^2+(y-sy12).^2 ));
      LLR = -(min( (x-sx0(i,:)).^2+(y-sy0(i,:)).^2 ) - min( (x-sx1(i,:)).^2+(y-sy1(i,:)).^2 ));
      L = [L LLR];
    end
    
%     L = [L LLR1 LLR2];
end  
%%%%%%

figure(1)
plot(normalize(L))
hold on
plot(normalize(demodSignal))

% figure(2)
% plot(demodSignal)
% plot(demodSignal,'.')
