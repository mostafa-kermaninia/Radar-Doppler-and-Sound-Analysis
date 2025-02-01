close all;
clc;
clear;
%% Part 1
fc = 5;
tstart = 0;
tend = 1;
fs = 100;

t = tstart: 1/fs : tend - 1/fs;
x = cos(2*pi*fc*t);
figure
plot(t, x)
xlabel('time');
ylabel('x(t)');
title('send signal');

%% Part 2

alpha = 0.5;
Beta = 0.3;
R = 250; %Km
V = 180; %Km / h
fd = Beta * V / 3.6;
c = 3*10^8;
td = 2/c * R * 1000;

y = alpha * cos(2*pi*(fc+fd)*(t-td));

figure
plot(t, y)
xlabel('time');
ylabel('y(t)');
title('recived signal');

%% Part 3
N = (tend - tstart) * fs;
f = -fs/2 : fs/N : fs/2-fs/N;

FTR = fftshift(fft(y));
FTR = FTR/max(abs(FTR));
[value idx1] = max(FTR(51:100));
freq = f(50+idx1);

phaseVal = abs(angle(FTR(50+idx1)));

fdnew = freq - fc;
tdnew = phaseVal/(2*pi*(fc+fdnew));
Vnew = fdnew * 3.6/Beta; 
Rnew = round(tdnew / 1000 * 0.5 * c);

disp(['fd = ', num2str(fdnew)]);
disp(['td = ', num2str(tdnew)]);
disp(['V = ', num2str(Vnew)]);
disp(['R = ', num2str(Rnew)]);

%% Part 4

std = 0.01;
noise = randn(1,length(y));   
True_R = 1;
True_V = 1;
for i=1:1000
    noise_i = std*noise;
    y_noisy = y + noise_i;
    
    FTR = fftshift(fft(y_noisy));
    FTR = FTR/max(abs(FTR));
    [value idx1] = max(FTR(51:100));
    freq = f(50+idx1);
    phaseVal = abs(angle(FTR(50+idx1)));

    fdnew = freq - fc;
    tdnew = phaseVal/(2*pi*(fc+fdnew));
    Vnew = fdnew * 3.6/Beta; 
    Rnew = round(tdnew / 1000 * 0.5 * c);
    
    if Rnew ~= R && True_R==1
        disp(['With std = ', num2str(std), 'Rnew = ', num2str(Rnew), ' detected wrong']);
        True_R=0;
    end
    
    if Vnew ~= V && True_V==1
        disp(['With std = ', num2str(std), 'Vnew = ', num2str(Vnew), ' detected wrong']);
        True_V=0;
    end
    
    if True_V + True_R == 0
        break;
    end
    std = std + 0.01;
end

%% Part 5

alpha1 = 0.5;
Beta = 0.3;
R1 = 250; %Km
V1 = 180; %Km / h
fd1 = Beta * V1 / 3.6;
td1 = 2/c * R1 * 1000;

alpha2 = 0.6;
R2 = 200; %Km
V2 = 216; %Km / h
fd2 = Beta * V2 / 3.6;
td2 = 2/c * R2 * 1000;

y1 = alpha1 * cos(2*pi*(fc+fd1)*(t-td1));
y2 = alpha2 * cos(2*pi*(fc+fd2)*(t-td2));
y = y1 + y2;

figure
plot(t, y)
xlabel('time');
ylabel('y(t)');
title('recived signal for 2 objects');

%% Part 6
FTR = fftshift(fft(y));
FTR = FTR/max(abs(FTR));

[value idx1] = max(FTR(51:100));
freq1 = f(50+idx1);
phaseVal1 = abs(angle(FTR(50+idx1)));
FTR(50+idx1) = 0;

[value idx2] = max(FTR(51:100));
freq2 = f(50+idx2);
phaseVal2 = abs(angle(FTR(50+idx2)));



fdnew1 = freq1 - fc;
tdnew1 = phaseVal1/(2*pi*(fc+fdnew1));
V1new = fdnew1 * 3.6/Beta; 
R1new = tdnew1 / 1000 * 0.5 * c;

fdnew2 = freq2 - fc;
tdnew2 = phaseVal2/(2*pi*(fc+fdnew2));
V2new = fdnew2 * 3.6/Beta; 
R2new = round(tdnew2 * c / 1000 * 0.5);

disp(['For object 1 V= ', num2str(V1new), 'and R= ', num2str(R1new)]);
disp(['For object 2 V= ', num2str(V2new), 'and R= ', num2str(R2new)]);


%% Part 7 test our Idea

alpha1 = 0.5;
Beta = 0.3;
R1 = 250; %Km
V1 = 180; %Km / h
fd1 = Beta * V1 / 3.6;
td1 = 2/c * R1 * 1000;

alpha2 = 0.6;
R2 = 200; %Km
V2 = 192; %Km / h
fd2 = Beta * V2 / 3.6;
td2 = 2/c * R2 * 1000;

y1 = alpha1 * cos(2*pi*(fc+fd1)*(t-td1));
y2 = alpha2 * cos(2*pi*(fc+fd2)*(t-td2));
y = y1 + y2;

FTR = fftshift(fft(y));
FTR = FTR/max(abs(FTR));

[value idx1] = max(FTR(51:100));
freq1 = f(50+idx1);
phaseVal1 = abs(angle(FTR(50+idx1)));
FTR(50+idx1) = 0;

[value idx2] = max(FTR(51:100));
freq2 = f(50+idx2);
phaseVal2 = abs(angle(FTR(50+idx2)));



fdnew1 = freq1 - fc;
tdnew1 = phaseVal1/(2*pi*(fc+fdnew1));
V1new = fdnew1 * 3.6/Beta; 
R1new = tdnew1 / 1000 * 0.5 * c;

fdnew2 = freq2 - fc;
tdnew2 = phaseVal2/(2*pi*(fc+fdnew2));
V2new = fdnew2 * 3.6/Beta; 
R2new = round(tdnew2 * c / 1000 * 0.5);

disp('test that if |V1-V2|>= 12 it works correctly');
disp(['For object 1 V= ', num2str(V1new), 'and R= ', num2str(R1new)]);
disp(['For object 2 V= ', num2str(V2new), 'and R= ', num2str(R2new)]);







