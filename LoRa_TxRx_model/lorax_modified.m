function [upchirp, check_chirp, downchirp, check_data] = lorax_modified( data, fs, B, Ts, SF, BW)
% �����������
N = length(data);       % ������ ������������������ ������

% ����� � �������
% fs = 1e6;              % ������� �������������
ts = 1/fs;              % ����� �������������
m = BW/Ts;              % ����������� �������� �������
t = 0:ts:Ts-ts;         % ������ ������� ������ �����
% t = 0:127;

% ������������ ������ ��������
code_word = zeros(1,SF); % ������� ������� �����
upchirp = [];           % �����

%% �������� ������
% upchirp = zeros(1,samples);
% chirp = A*cos(2*pi*(0*t + (m*t.^2)/2));
% downchirp = flip(chirp);
% check_chirp = (1/sqrt(B))*exp( 1i * (2*pi*(0*t + (m*t.^2)/2)) );
% chirp = (1/sqrt(B))*exp( 1i * (2*pi*(0*t + (m*t.^2)/2)) );
% downchirp = (1/sqrt(B))*exp( -1i * (2*pi*(0*t + (m*t.^2)/2)) );

check_chirp = exp( 1i * (2*pi*(0*t + (m*t.^2)/2)) );
chirp = exp( 1i * (2*pi*(0*t + (m*t.^2)/2)) );
downchirp = exp( -1i * (2*pi*(0*t + (m*t.^2)/2)) );

%% ���������
for i = 1:N/SF
code_word = bi2de(data(SF*i-SF+1:SF*i));      % �������� �������� �����
check_data(i) = code_word; % �������� ��������������� �������
cs = single((code_word/B)*Ts/ts);             % ����� ������

chirp1 = chirp(cs+1:end); % ������ ����
chirp2 = chirp(1:cs); % ������ ����

%%%%%%%
upchirp = [upchirp, chirp1, chirp2];          % �������� ����
end

end