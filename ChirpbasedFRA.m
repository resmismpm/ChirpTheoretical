function z_singlechirp = ChirpbasedFRA(sys,Index)

syms s t1

% Convert transfer function to symbolic function
[Num,Den] = tfdata(sys);
G = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s);

%% Input Signal
tfinal = 10;
t = linspace(0,tfinal,250000);

if Index == 1  %% Linear Chirp    
    f0 = 1 ; f1 = 400 ; Phi0 = 0;
    c = (f1-f0)/tfinal;
    Phi = Phi0+ 2*pi*(0.5*c*t1.^2 + f0*t1);
    
elseif Index == 2  %% Exponential Chirp
    f0 = 1; f1 = 1000 ;
    Phi0 = 0;
    c = (f1/f0)^(1/tfinal);
    Phi = Phi0+ 2*pi*f0*(c.^t1 - 1)/log(c);
    
elseif Index == 3 %% 4th order Chirp
    f0 = 1 ; f1 = 1000;
    Phi0 = 0;
    c = (f1-f0)/(5*tfinal^4);
    Phi = Phi0+ 2*pi*(c*t1.^5 + f0*t1);
end

Omega = diff(Phi);
Omega_True = matlabFunction(Omega);
Omega_TrueData = Omega_True(t);
Phi_True = matlabFunction(Phi);
Phi_TrueData = Phi_True(t);

%% Theoretical frequency response of the system

AR_True = matlabFunction(abs(subs(G,s,1i*Omega)));
AR_TrueData = AR_True(t);
PhaseLag_True = matlabFunction(phase(subs(G,s,1i*Omega)));
PhaseLag_TrueData = PhaseLag_True(t);
z_True = AR_TrueData.*exp(1i*PhaseLag_TrueData);
x_pred = AR_TrueData.*sin(Phi_TrueData+PhaseLag_TrueData);
x_pred = x_pred';

%% Simulation of transfer function
u = sin(Phi_TrueData);
x=lsim(sys,u,t);

figure
plot(t,x,'LineWidth',2,'color',[0.47, 0.67, 0.19])
hold all
plot(t,x_pred,'-.','LineWidth',1,'color',[0, 0.45, 0.74])
plot(t,AR_TrueData,'--','LineWidth',1,'color',[0.64, 0.08, 0.18])
plot(t,-AR_TrueData,'--','LineWidth',1,'color',[0.64, 0.08, 0.18])
% xlim([0 0.6])
xlabel('Time (sec)')
ylabel('x')
legend('True chirp response','Assumed response','True AR','FontSize',10,'Orientation','horizontal')
title(['G(s) =', arrayfun(@char, G, 'uniform', 0)],'FontSize',12)

saveas(gcf,'Figure1.png')


% figure
% plot(t,(x-x_pred)./x,'LineWidth',2)
% hold all
% title(['G(s) =', arrayfun(@char, G, 'uniform', 0)],'FontSize',12)
% xlabel('Time (sec)')
% ylabel('Error %')
figure
plot(t,(x-x_pred),'LineWidth',2)
hold all
title(['G(s) =', arrayfun(@char, G, 'uniform', 0)],'FontSize',12)
xlabel('Time (sec)')
ylabel('Error')

saveas(gcf,'Figure2.png')

%% Chirp-based FRA
if Index ==1
    no_s = round(0.03*length(t));
else
    no_s = round(0.12*length(t));
end
t = t(no_s:end);
x    = x(no_s:end);

Omega_TrueData = Omega_TrueData(no_s:end);
Phi_TrueData = Phi_TrueData(no_s:end);
AR_TrueData = AR_TrueData(no_s:end);
PhaseLag_TrueData = PhaseLag_TrueData(no_s:end);
z_True = z_True(no_s:end);

%% Finding amplitude ratio using envelope of output chirp signal
[A_singlechirp,~] = envelope(x,round(6e-4*length(t)),'peak');

%% Finding phase lag
A_singlechirp1 = A_singlechirp*max(abs(x./A_singlechirp));  % To ensure none of the values are greater than one
x_scaled = (x./A_singlechirp1);
phi_half1 = acos(-x_scaled);

phi_half3 = phi_half1;
t2 = t;

complex_indices = find(imag(phi_half3)~=0);
phi_half3(complex_indices) = [];
t2(complex_indices) = [];
phi_half1(complex_indices) = spline(t2,phi_half3,complex_indices);

% Note that here acos(-x) is used even though the input is sin signal. But this
% difference is accounted by subtracting a pi/2 while calculating the final
% phase lag. This is coz sin(Phi) = -cos(Phi-pi/2)

%% Converting phi_half1 in [0,pi] to [-pi, pi] & then to a continuously increasing phase using unwrap

[max_phi(:,2), max_phi(:,1)]=findpeaks(phi_half1);                    %% Finding peaks of v1
[min_phi(:,2), min_phi(:,1)]=findpeaks(-phi_half1);                    %% Finding peaks of v1

phi_half2 = phi_half1;

if max_phi(1,1) < min_phi(1,1)
    for jdx = 1: min(size(max_phi,1),size(min_phi,1))
        phi_half2(max_phi(jdx,1):min_phi(jdx,1)) = -phi_half2(max_phi(jdx,1):min_phi(jdx,1));
    end
else
    for jdx = 1: min(size(max_phi,1),size(min_phi,1))-1
        phi_half2(max_phi(jdx,1):min_phi(jdx+1,1)) = -phi_half2(max_phi(jdx,1):min_phi(jdx+1,1));
    end
end
phi_half3_snglechrp = unwrap(phi_half2);

%% Calculating phase lag and impedance
phi_singlechirp = phi_half3_snglechrp-Phi_TrueData'-pi/2;
phi_singlechirp = csaps(t,phi_singlechirp,0.8,t);
z_singlechirp = A_singlechirp .* exp(phi_singlechirp' *1i);

%% Plots
figure
plot(real(z_singlechirp),imag(z_singlechirp),'color',[0.47, 0.67, 0.19]);
hold all
plot(real(z_True),imag(z_True),'-.','color',[0.64, 0.08, 0.18]);
legend('Chirp Analysis','Theoretical','Orientation','Horizontal')
xlabel('Re(z)');ylabel('Im(z)')
title(['G(s) =', arrayfun(@char, G, 'uniform', 0)],'FontSize',12)
saveas(gcf,'Figure3.png')

end
