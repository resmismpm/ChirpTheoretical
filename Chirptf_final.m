function z_singlechirp = Chirptf_final(sys,tfinal, f0, f1,ts)
% f0 = 1e-3;  f1 = 1e2;
% ts = 1e-3;
% 
% tfinal = 200;
t = 0:ts:tfinal;

a = (f1-f0)/(5*tfinal^4);
Phi = 2*pi*(a*t.^5 + f0*t);
f  = 5*a*t.^4 + f0;
% c = (f1/f0)^(1/tfinal);
% Phi = 2*pi*f0*(c.^t - 1)/log(c);
x = sin(Phi);

y = lsim(sys,x,t);


%% Single Chirp Analysis

%% Finding amplitude ratio using envelope of output chirp signal

[maxtab1, ~] = peakdet(y,0.0002);
if maxtab1(1,1) ~=1
    maxtab1 = [1 maxtab1(1,2); maxtab1];
end
A_singlechirp = spline(maxtab1(:,1),maxtab1(:,2),1:length(y))';  %% fitting a spline to get max values at each time instant
% A_singlechirp = smooth(A_singlechirp,length(y)/50,'moving');

% [A_singlechirp,~] = envelope(y,1000,'peak');

%% Removing the transients
no_s = round(0.05*length(t));
no_f = 0;%round(0.001*length(t));
t = t(no_s:end-no_f);
f = f(no_s:end-no_f);
x = x(no_s:end-no_f);
y = y(no_s:end-no_f);
Phi = Phi(no_s:end-no_f);
A_singlechirp = A_singlechirp(no_s:end-no_f);

%% Finding phase lag
% A_singlechirp1 = A_singlechirp*max(abs(y./A_singlechirp));  % To ensure none of the values are greater than one
y_scaled = (y./A_singlechirp);
phi_half1 = acos(-y_scaled);
ind = find(imag(phi_half1)~=0);
phi_half3 = phi_half1;
phi_half3(ind) = [];
xind = 1:length(phi_half1);
xind(ind) = [];
phi_half1(ind) = spline(xind,phi_half3,ind);

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
phi_singlechirp = phi_half3_snglechrp-Phi'-pi/2;
phi_singlechirp = csaps(t,phi_singlechirp,0.6,t);
z_singlechirp = A_singlechirp .* exp(phi_singlechirp' *1i);
plot(z_singlechirp,'.')
% hold all

end
% x_fft = fft(x);
% N = length(x);
% x_fft = x_fft(1:end/2);
% 
% F = ((1:N/2)/N)*(1/ts);

% stem(F,x_fft/max(x_fft))









