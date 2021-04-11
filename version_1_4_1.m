clear all;
cla('reset');
close all;

NmberOfPeriods=10;% ���-�� ��������
amount_elem=2;
T=100*10^-3;% (T=10^-1;�� �������) % ������ �������
C=3*10^8; % �������� ����� � �����
f0=24.15*10^9; % ������� �������
dF=200*10^6; %2*�������� �������
Max_dist=150;
range_resolution=1.5;

Ampl=1; % ���������
sigma=0.1*Ampl;% ��������� ����
step=10^3;% kol-vo tochek na periud
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Object=3;

auto_x=5;%  ��������� ����������
auto_y=0;% ��������� ����������

speed_bike_x=0; % �������� ������������� �� ��� � ������ � ������� = 15 ��/�
speed_bike_y=0;% �������� ������������� �� ��� y ������ � ������� = 1 ��/�
speed_auto_x=20; % �������� ���� ������ �� ��� � �/� = 60 ��/�
speed_auto_y=0; % �������� ���� ������ �� ��� y �/� = 60 ��/�

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt=T/step;% ��������� ���
N=NmberOfPeriods*step;%kol-vo tochek
speed_dif_x=speed_auto_x-speed_bike_x; %���������� �������� �� ��� x
speed_dif_y=speed_auto_y-speed_bike_y; %���������� �������� �� ��� y

% �������� ������������� �������
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
% title('����� ��������� ������� ����������� �������');


time_dalay=zeros(1,N);% �������� ���������� ������� � ������
angle_auto=zeros(1,N);% ���� �� ������� ��������� ����
dif_fasa=zeros(1,N);% �������� ��� ����� ��������� ���������� �������� ���
n=1;
while n<=N
    radius_to_auto=sqrt(auto_x^2+auto_y^2);% ������ �� �������
    time_dalay(n)=(2*radius_to_auto)/C;% ����� �������� �������
    angle_auto=atan(auto_y/auto_x)*180/pi;% � ��������
    dif_fasa(n)=2*pi*0.5*sin(angle_auto*pi/180);
    auto_x=auto_x+dt*speed_dif_x; % ��������� ��������� ���� �� ��� �
    auto_y=auto_y+dt*speed_dif_y;% ��������� ��������� ���� �� ��� y
    n=n+1; 
end




% ��������� ��������� ������

saw_in=zeros(1,N);% ��������� ��������� ������
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
 title('������� ����������� � ��������� �������');


%������ �������� ����������� � ���������� ��������
saw_dif=zeros(1,N);
n=1;
while n<N
    saw_dif(n)=saw_in(n)-saw_out(n);
    n=n+1;
end

figure;
plot(abs(saw_dif),'b');
title('���������� �������');



dif_freq_theory=(2*radius_to_auto*dF*2)/(C*T); % ���������� ������� ��c��������� ������������
dif_freq_resalt=max(saw_dif);% ���������� ������� ��c��������� �� ��������

% �������� �������� �� ���������� ������� � �����
speed_auto_mod=zeros(1,N);% ����������������� ��������� �� ����
radius_to_auto_mod=zeros(1,N);


fsamp=2*2667;% ������� ������������� ��������
d=1/(fsamp); %��� ��������������

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

 %�������� ���������� ������� �� ���������� �������
 nomber=char(48+item_number);
 Name_faile=strcat(nomber,'_elem','.wav');
 load handel.mat;
 Fs=round(fsamp);
 audiowrite(Name_faile,signal_dif_freq,Fs);
 info = audioinfo(Name_faile);
end

figure;
plot(signal_dif_freq);
title('������������� ������ �� ���������� �������');

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
title('������ �������');
i=i+part;
end

