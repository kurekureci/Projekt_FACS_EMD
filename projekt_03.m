close all; clear all; clc;

%% nacteni a uprava vstupniho signalu (pro porovnani RMSE)
load('202m.mat');
x = val(1,:)-mean(val(1,:));
d = length(x);

N = 400;
fvz = 360;  % z databaze MIT-BIH Arrhythmia
Wnv = ([0.7/(fvz/2),100/(fvz/2)]); % filtrace driftu
bv = fir1(N,Wnv,'bandpass'); % výpoèet imp. char. filtru
yv = filtfilt(bv,1,x); % vysledek filtrace
x = yv;
Wnv2 = ([55/(fvz/2),65/(fvz/2)]); % pro filtraci brumu 60 Hz
bv2 = fir1(N,Wnv2,'stop'); % výpoèet imp. char. filtru
yv2 = filtfilt(bv2,1,x); % vysledek filtrace
x = yv2;

%% pridani brumu a driftu
A = 50; % amplituda brumu (uV)
f = 50; % frekvence brumu(Hz)
cas = (1:d)/fvz; 
brum = A*sin(2*pi*f*cas); % sitove ruseni
xBrum = x + brum;

A2 = 150; % amplituda driftu (uV)
f2 = 0.18; % frekvence driftu (Hz)
drift = A2*sin(2*pi*f2*cas);
xDrift = x + drift;
xSum = xBrum + drift;

figure(1);
subplot(3,1,1)
plot(x,'linewidth',2);
title('Vstupní EKG signál','FontSize',17);
xlim([0 d]);
subplot(3,1,2)
plot(xBrum,'linewidth',2);
title(['Pøidání brumu (frekvence: 50 Hz, amplituda: ' num2str(A) ' \muV)'],'FontSize',17);
ylabel('napìtí [\muV]','FontSize',15)
xlim([0 d]);
subplot(3,1,3)
plot(xSum,'linewidth',2);
title(['Pøidání driftu (frekvence: ' num2str(f2) ' Hz, amplituda ' num2str(A2) ' \muV)'],'FontSize',17);
xlabel('poøadí vzorku [-]','FontSize',15)
xlim([0 d]);

%% EMD
Y = xSum; % vstupni signal
NoiseLevel = 0; % uroven sumu
NE = 1; % pro EMD
numImf = -1; % pocet imf
imfs = eemd(Y, NoiseLevel, NE, numImf);

%% zobrazeni IMF
[a,b]=size(imfs);
pul = ceil(a/2);

figure(2);
for i = 1:pul
   subplot (pul,1,i);
   plot(imfs(i,:),'linewidth',2);
   ylabel(['IMF ' num2str(i)],'FontSize',15);
   xlim([0 b]);
   
   if i == 1
       title('Výsledek EMD, vnitøní funkce signálu','FontSize',17);
   end
   
end
xlabel('poøadí vzorku [-]','FontSize',15)
c = 1; 

figure(3);
for i = pul+1:a
   subplot (floor(a/2),1,c);
   plot(imfs(i,:),'linewidth',2);
   ylabel(['IMF ' num2str(i)],'FontSize',15);
   c = c+1;
   xlim([0 b]);
   
   if i == pul+1
       title('Výsledek EMD, vnitøní funkce signálu','FontSize',17);
   end
   
end
xlabel('poøadí vzorku [-]','FontSize',15)

%% Odstraneni sitoveho ruseni - filtrace 1. IMF FIR filtrem DP s mezni frekvenci 31 Hz 
N = 400; % délka impulzní charakteristiky filtru
Wn = (31/(fvz/2)); % mez
b1 = fir1(N,Wn,'low'); % výpoèet imp. char. filtru
y1 = filtfilt(b1,1,imfs(1,:)); % vysledek filtrace
bezBrumu = xSum - imfs(1,:) + y1;

figure(4);
subplot(3,1,1);
plot(imfs(1,:),'linewidth',2);
title('IMF 1','FontSize',17);
xlim([0 d]);
subplot(3,1,2);
plot(y1,'linewidth',2);
title('1. IMF po filtraci dolní propustí s mezní frekvencí 31 Hz','FontSize',17);
ylabel('napìtí [\muV]','FontSize',15)
xlim([0 d]);
subplot(3,1,3);
plot(bezBrumu,'linewidth',2);
title('Výsledek odstranìní brumu','FontSize',17);
xlabel('poøadí vzorku [-]','FontSize',15)
xlim([0 d]);

%% Odstraneni driftu - zjisteni poctu ovlivnenych IMF
B = zeros(1,b); % b = pocet vzorku

for k = a:-1:1 
    spektrum = abs(fft(imfs(k,:))); % vypocet amplitudoveho spektra IMF
    [m,ind] = max(spektrum);
    frekv = ind(end) * (fvz/b); % zisk frekvence IMF
    
    if frekv < 0.5
        B = B + imfs(k,:);
    end
end

bezDriftu = xSum - B;
bezSumu = bezBrumu - B;

figure(5);
subplot(3,1,1);
plot(bezBrumu,'Color','b','linewidth',2);
title('EKG signál s driftem','FontSize',17);
xlim([0 d]);
subplot(3,1,2);
plot(B,'linewidth',2);
title('Souèet IMF ovlivnìných driftem','FontSize',17);
ylabel('napìtí [\muV]','FontSize',15)
xlim([0 d]);
subplot(3,1,3);
plot(bezSumu,'linewidth',2);
title('Výsledek po odstranìní driftu','FontSize',17)
xlabel('poøadí vzorku [-]','FontSize',15)
xlim([0 d]);

%% vypocet stredni kvadraticke odchylky RMSE 
mse = (sum((x-bezSumu).^2))/b;
rmse = sqrt(mse);
disp(['RMSE: ' num2str(rmse)]);

%% vypocet SNR
w = brum + drift;
snrIn = 10*log10((sum(x.^2))/(sum(w.^2)));
snrOut = 10*log10((sum(x.^2))/(sum((bezSumu-x).^2)));
disp(['SNRin: ' num2str(snrIn) ' dB']);
disp(['SNRout: ' num2str(snrOut) ' dB']);

%% Srovnani vysledku s originalem
figure(6);
plot(x,'Color','b','linewidth',2);
hold on;
plot(bezSumu,'Color','r','linewidth',2);
title('Porovnání výsledku filtrace se vstupním signálem EKG' ,'FontSize',17);
leg = legend('Vstupní EKG','Výsledné EKG');
set(leg,'Location','NorthEast')
set(leg,'FontSize',12);
ylabel('napìtí [\muV]','FontSize',15)
xlabel('poøadí vzorku [-]','FontSize',15)
xlim([0 d]);
