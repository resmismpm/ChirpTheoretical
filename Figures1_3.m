clc
clear

%% Define the transfer function of interest
sys = tf(2,[0.01 1]);
% sys = tf(40,[1 400 40000]);
% sys = tf(2,[1/10000 1/50 1]);
% sys = tf(2,[1 1000 420000 112000000 16000000000]);
% sys = tf([2 40],[1 130 3000]);
% sys = tf(2,[1 400 380000  68000000 28900000000]);


%% Chirp-based FRA

ChirpIndex = 3; % 1 for Linear, 2 for Exponential and 3 for Fourth-order input chirp signal

z_singlechirp = ChirpbasedFRA(sys,ChirpIndex);