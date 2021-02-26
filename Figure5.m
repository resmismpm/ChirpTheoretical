%%% Code used for generating the experimental result in Nature
%%% Computational Science paper draft
%%% Final code for dual chirp calculation for hardware data
%%% Comparison of hardware results with theoretical Nyquist plot

%% Hardware Calculations

%% Step 1 :  Load data for the corresponding circuit

clc
clear

load('HardwareData.mat')

%% Step 2: Hardware data postprossessing

jdx=5000;
while jdx < length(V1out)
    if abs(V1out(jdx) - V1out(jdx-1))>0.15
        cut_1 = jdx-1;
        break
    else
        jdx = jdx+1;
    end
end

idx=5000;
while idx < length(V2out)
    if abs(V2out(idx) - V2out(idx-1))>0.1
        cut_2 = idx-1;
        break
    else
        idx = idx+1;
    end
end

Vout_1 = V1out(cut_1 + 1:80000+cut_1);
Vout_2 = V2out(cut_2 + 1:80000+cut_2);

v1_2 = Vout_1;
v2_2 = Vout_2;
tout = t;

%% Step 2 (a): Denoising the voltage signals

%% (i) Using inbuilt wavelet denoising upto 13000 data points
[v1_den1,cfs11] = cmddenoise(v1_2(1:42000,1),'sym4',24,'s');
v1_den1 = smooth(v1_den1,80,'moving');
[v1_den2,cfs12] = cmddenoise(v1_2(40001:65000,1),'sym4',2,'s');
v1_den = [v1_den1(1:41000); v1_den2(1001:end)'; v1_2(65001:end)];

[v2_den1,cfs21] = cmddenoise(v2_2(1:42000,1),'sym4',24,'s');
v2_den1 = smooth(v2_den1,80,'moving');
[v2_den2,cfs22] = cmddenoise(v2_2(40001:65000,1),'sym4',2,'s');
v2_den = [v2_den1(1:41000); v2_den2(1001:end)'; v2_2(65001:end)];

v1_2 = v1_den;
v2_2 = v2_den;

%% Step 3 : Dual chirp calculation

%% Step 3 (a): Frequecy calculation using exponential chirp

f0 = 0.1;
ff = 10000;
tf = 2;

freq = f0*((ff/f0).^(tout'/tf));
phi_rpi = 3*pi/2 + 2*pi*f0*((freq/f0)-1)./(log((ff/f0).^(1/tf)));

%% Step 3 (b): Removing few data points from start and end

no_s = 2500;        % No. of points deleted from starting
no_f = 60;        % No. of points deleted from end

v1_2 = v1_2(no_s : end-no_f);
v2_2 = v2_2(no_s : end-no_f);
phi_rpi = phi_rpi(no_s : end-no_f);
freq = freq(no_s : end-no_f);
tout = tout(no_s : end-no_f);
%% Step 3 (c): Amplitude and phase calculation using dual chirp equation

phi_half1_rpi = (atan((cos(pi/6)-v2_2./v1_2)/sin(pi/6)));
[maxphi_rpi, minphi_rpi] = peakdet(phi_half1_rpi,1);

phi_half_rpi = phi_half1_rpi;

j=1;    idx = 0;
ind_old =1;
while j <= length(maxphi_rpi)
    ind = maxphi_rpi(j,1);
    phi_half_rpi(ind_old:ind) = phi_half1_rpi(ind_old:ind) + ((j+idx-1)*pi);
    ind_old = ind+1;
    
    if j <= length(minphi_rpi)
        ind = minphi_rpi(j,1);
        if ind ~= ind_old
            for m = ind_old : ind-1
                idx = idx+1;
                phi_half_rpi(m) = phi_half1_rpi(m) + ((j+idx-1)*pi);
            end
            ind_old = ind;
        end
    end
    j = j+1;
end
phi_half_rpi(ind_old:end) = phi_half1_rpi(ind_old:end) + (j+idx-1)*pi ;

psi1_rpi = phi_half_rpi - phi_rpi;

A_rpi = abs(v1_2./cos(phi_rpi + psi1_rpi));


%% Step 3 (4): Smoothing of amplitude and psi

% Phase
psi1_rpi_fit =smooth(psi1_rpi,5000,'moving');

% Amplitude
A_rpi_fit =smooth(A_rpi,5000,'moving');

%% Step 3 (d): Impedance from hardware data

z_rpi =  R_series./(((exp(-1i*psi1_rpi_fit))./A_rpi_fit)-1);



%% Step 4: Comparison of the impedance from matlab and Rpi

figure
plot(real(z_rpi),imag(z_rpi),'.')
hold all
plot(ZData(:,2))

xlabel('Re(Z)')
ylabel('Im(Z)')
title('Impedance Plot')
legend('Dual chirp - experimental' ,'Theoretical')
saveas(gcf,'Figure5.png')

