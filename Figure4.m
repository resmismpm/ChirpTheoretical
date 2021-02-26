clc
clear
close all

%% System of interest
sys = tf([2 40],[1 130 300]);

%% Interested frequency regime and signal duration
tfinal = 1000;   % Signal duration

f0 = 1e-3;   % Initial Frequency
f1 = 1e2;    % Final Frequency
ts = 1e-3;   % Sampling Time

%% Theoretical FRA
[H,~] = freqresp(sys,logspace(log10(f0),log10(f1),1000),'Hz');
z_true = squeeze(H);
plot(z_true,'--')
hold all

%% Chirp-based FRA
z_singlechirp = Chirptf_final(sys,tfinal, f0, f1,ts);

%% Sum-of-sines FRA
n = 1000;
z = sumofsines(sys,tfinal, f0, f1,ts,n);

legend('Theoretical Nyquist plot','Chirp signal-based FRA', 'Multi-sine FRA')
saveas(gcf,'Figure4.png')
