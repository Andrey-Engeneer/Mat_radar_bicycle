clear all;
cla('reset');
close all;

NmberOfPeriods=10;% кол-во периудов
amount_elem=2;
T=100*10^-3;% (T=10^-1;по заданию) % периуд сигнала
C=3*10^8; % скорость света в среде
f0=24.15*10^9; % несущая частота
dF=200*10^6; %2*дивиация частоты
Max_dist=150;
range_resolution=1.5;

Ampl=1; % амплитуда
sigma=0.1*Ampl;% дисперсия шума
step=10^3;% kol-vo tochek na periud
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Object=3;

auto_x=5;%  начальные координаты
auto_y=0;% начальные координаты

speed_bike_x=0; % скорость велосипедиста по оси х метров в секунду = 15 км/ч
speed_bike_y=0;% скорость велосипедиста по оси y метров в секунду = 1 км/ч
speed_auto_x=20; % скорость авто метров по оси х м/с = 60 км/ч
speed_auto_y=0; % скорость авто метров по оси y м/с = 60 км/ч

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt=T/step;% временной шаг
N=NmberOfPeriods*step;%kol-vo tochek
speed_dif_x=speed_auto_x-speed_bike_x; %разностная скорость по оси x
speed_dif_y=speed_auto_y-speed_bike_y; %разностная скорость по оси y

% создание пилообразного сигнала
saw_out=zeros(1,N);
i=0;

while i<N
n=1;
    while n*dt<=T
        if n*dt<T/2
            saw_out(i+n)=(dF)/(T/2)*n*dt+(f0-dF/2);
                else
            saw_out(i+n)=(dF)/(T/2)*(T-n*dt)+(f0-dF/2);
        end
        n=n+1;
    end
    i=i+n-1; 
end 

% figure;
% plot(saw_out,'r');
% title('закон изменения частоты излучённого сигнала');


time_dalay=zeros(1,N);% задержка отражённого сигнала в тактах
angle_auto=zeros(1,N);% угол на котором находится авто
dif_fasa=zeros(1,N);% разность фаз между соседними элементами антенной реш
n=1;
while n<=N
    radius_to_auto=sqrt(auto_x^2+auto_y^2);% радиус до объекта
    time_dalay(n)=(2*radius_to_auto)/C;% время задержки сигнала
    angle_auto=atan(auto_y/auto_x)*180/pi;% в градусах
    dif_fasa(n)=2*pi*0.5*sin(angle_auto*pi/180);
    auto_x=auto_x+dt*speed_dif_x; % изменение координат авто по оси х
    auto_y=auto_y+dt*speed_dif_y;% изменение координат авто по оси y
    n=n+1; 
end




% пришедший отражённый сигнал

saw_in=zeros(1,N);% пришедший отражённый сигнал
i=0;
while i<N
n=1;
    while n*dt<=T
        if (n*dt)<time_dalay(i+n)
           saw_in(i+n)=(dF)/(T/2)*(time_dalay(i+n)-n*dt)+(f0-dF/2)+2*speed_dif_x/(saw_out(n)/C); 
        elseif(n*dt)<=T/2+time_dalay(i+n)
            saw_in(i+n)=(dF)/(T/2)*(n*dt-time_dalay(i+n))+(f0-dF/2)+2*speed_dif_x/(saw_out(n)/C);
        elseif(n*dt)>T/2+time_dalay(i+n)
            saw_in(i+n)=(dF)/(T/2)*(T+time_dalay(i+n)-n*dt)+(f0-dF/2)+2*speed_dif_x/(saw_out(n)/C);
        end
        n=n+1;
    end
    i=i+n-1; 
end 
 figure;
 plot(saw_in,'b');
 hold on;
 plot(saw_out,'r');
 hold off;
% 
 title('Частота излучённого и принятого сигнала');


%график разности излучённого и отражённого сигналов
saw_dif=zeros(1,N);
n=1;
while n<N
    saw_dif(n)=saw_in(n)-saw_out(n);
    n=n+1;
end

figure;
plot(abs(saw_dif),'b');
title('Разностная частота');



dif_freq_theory=(2*radius_to_auto*dF*2)/(C*T); % разностная частота раcсчитанная теоретически
dif_freq_resalt=max(saw_dif);% разностная частота раcсчитанная из графиков

% создание косинуса на разностной частоте с шумом
speed_auto_mod=zeros(1,N);% промоделированное растояние до авто
radius_to_auto_mod=zeros(1,N);


fsamp=2*2667;% частота дискретизации косинуса
d=1/(fsamp); %шаг дискретизацции

lasting_signal=N*dt;
quantity_dots=round(lasting_signal/d);
signal_dif_freq=zeros(1,quantity_dots);
noise_dif_freq=randn(1,quantity_dots);

test=zeros(1,quantity_dots);

for item_number=0:amount_elem-1
i=1;
n=1;
while n<quantity_dots
    if i*dt<=d*n
	i=i+1;
    end
    signal_dif_freq(n)=cosd(n*d*360*(abs(saw_dif(i))+dif_fasa(i)*item_number));
    signal_dif_freq(n)=Ampl*signal_dif_freq(n)+sigma*noise_dif_freq(n);
   
    n=n+1;
end

 %создание аудиофайла сигнала на разностной частоте
 nomber=char(48+item_number);
 Name_faile=strcat(nomber,'_elem','.wav');
 load handel.mat;
 Fs=round(fsamp);
 audiowrite(Name_faile,signal_dif_freq,Fs);
 info = audioinfo(Name_faile);
end

figure;
plot(signal_dif_freq);
title('Гармонический сигнал на разностной частоте');

[y,FFF] = audioread('0_elem.wav');

part=256;
i=0;
[m,n]=size(y);
part_signal=zeros(1,part);
while i<1%m 
    for n=1:part
    part_signal(n)=y(n+i);
    end
FftL=part;
spectr=abs(fft(part_signal,FftL));%*2/FftL;
F=0:(1/d)/FftL:(1/d)/2-1/FftL;
figure;
plot(F,abs(spectr(1:length(F))));
title('спектр сигнала');
i=i+part;
end

