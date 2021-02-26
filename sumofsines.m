function z = sumofsines(sys,tfinal, f0, f1,ts,n)

f = logspace(log10(f0),log10(f1),n);

%% SumofSines

t = 0:ts:tfinal;

x = zeros(n,length(t));

for k = 1: n
    x(k,:) = sin(2*pi*f(k)*t);
end

Input = sum(x);

%% Simulation

y = lsim(sys,Input,t);

N = length(Input);
x_fft = fft(Input');
x_fft = x_fft(1:length(x_fft)/2);

F = ((1:N/2)/N)*(1/ts);
[~, ind] = sort(x_fft,'descend');
xfftnew = x_fft(ind(1:n));

y_fft = fft(y);
y_fft = y_fft(1:length(x_fft)/2);
yfftnew = y_fft(ind(1:n));
z = yfftnew./xfftnew;

plot(z,'o')
end
