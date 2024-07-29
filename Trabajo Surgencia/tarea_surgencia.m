clc
clear all

load("VientosCosteros.mat");
%V_i componente Norte
%U_i componente Este

%Velocidad de los vientos en 37°S, 30°S y 21°S

rho_aire = 1.2; % kg/m^3
C_d = 1.3e-3;
rho_agua = 1025; % kg/m^3
omega = 7.292e-5; % 1/s
lat = [37 30 21];
costa = [27 10 4];

for i = 1:3
    angulo(:,i) = atand(Ui(:,i)./Vi(:,i)); %Angulo c/r al norte (eje y)
    
            %Debemos hayar los vectores de viento cuyas componentes apunten hacia el
            %"sur", de modo que
    
    A = find(Ui(:,i) > 0 & Vi(:,i) < 0);
    B = find(Ui(:,i) <= 0 & Vi(:,i) < 0);
    angulo(A,i) = angulo(A,i) + 180; %Rotar en angulo para que tenga componente norte
    angulo(B,i) = angulo(B,i) - 180;
    
            %Ahora consideramos el angulo que tiene la costa en la latitud dependiendo
            %con respecto al norte. Este valor debe restarse al angulo medido para
            %considerar adecuadamente la costa
    
    angulo(:,i) = angulo(:,i) - costa(i);
    clear A B
            %determinamos la magnitud del vector
    
    magnitud(:,i) = sqrt((Vi(:,i).^2) + (Ui(:,i).^2));
    f(i) = 2*omega*sind(-lat(i));
    for j = 1:3769
        y(j,i) = magnitud(j,i) * cosd(angulo(j,i));
        x(j,i) = magnitud(j,i) * sind(angulo(j,i));
        tau_y(j,i) = rho_aire * C_d * y(j,i) * magnitud(j,i); % Pa
        M_x(j,i) = tau_y(j,i) / (rho_agua*f(i)); % m^3/s/m
    end
end
%% 
viento = [2 5 10 15 20 25]; % m/s

CD_aleatorios = linspace(0.9e-3, 2.5e-3, 6);

for i = 1:6
    for j = 1:6
        tau(i,j) = rho_aire * CD_aleatorios(i) * viento(j).^2; % Pa
    end
end

%primera columna es para 5 m/s y seis distintos Cd

figure()
hold on
for i = 1:6
    plot(CD_aleatorios',tau(:,i),'LineWidth',2)
end
axis tight
grid minor
xlabel("Coeficiente de arraste C_D")
ylabel("Esfuerzo del viento \tau [Pa]")
title("Coeficientes de arrastre aleatorios y velocidades de viento.")
legend(string(viento + " m/s"),'Location','northwest')

%%
fecha = datetime(fecha);
figure()
    plot(fecha,M_x(:,1),'b','LineWidth',2)
    hold on
    yline(0,'LineWidth',2)
    xlabel("Fecha",'FontSize',20)
    ylabel("Transporte de Ekman [m^3/s/m]",'FontSize',20)
    title("Transporte de Ekman zonal 37°S",'FontSize',25,'FontWeight','bold')
    grid minor
figure()
    plot(fecha,M_x(:,2),'r','LineWidth',2)
    hold on
    yline(0,'LineWidth',2)
    xlabel("Fecha",'FontSize',20)
    ylabel("Transporte de Ekman [m^3/s/m]",'FontSize',20)
    title("Transporte de Ekman zonal 30°S",'FontSize',25,'FontWeight','bold')
    grid minor
figure()
    plot(fecha,M_x(:,3),'Color',[0.9290 0.6940 0.1250],'LineWidth',2)
    hold on
    yline(0,'LineWidth',2)
    xlabel("Fecha",'FontSize',20)
    ylabel("Transporte de Ekman [m^3/s/m]",'FontSize',20)
    title("Transporte de Ekman zonal 21°S",'FontSize',25,'FontWeight','bold')
    grid minor


%%
figure()
subplot(1,3,1)
    histogram(M_x(:,3),'FaceColor',[0.9290 0.6940 0.1250])
    title("21°S",'FontSize',15,'FontWeight','bold')
    xlabel("Transporte de Ekman [m^3/s/m]")
    annotation('textbox',[.13 .77 .1 .2], 'String','a)','EdgeColor','none')
    grid minor
subplot(1,3,2)
    histogram(M_x(:,2),'FaceColor','r')
    title("30°S",'FontSize',15,'FontWeight','bold')
    xlabel("Transporte de Ekman [m^3/s/m]")
    annotation('textbox',[.41 .77 .1 .2], 'String','b)','EdgeColor','none')
    grid minor
subplot(1,3,3)
    histogram(M_x(:,1),'FaceColor','b')
    title("37°S",'FontSize',15,'FontWeight','bold')
    xlabel("Transporte de Ekman [m^3/s/m]")
    annotation('textbox',[.7 .77 .1 .2], 'String','c)','EdgeColor','none')
    grid minor
    
%% 
fecha = datevec(fecha);
count = 0;
for i = 1:12
    A = find(fecha(:,2)==i);
    count = count+1;
    for j = 1:3
        dat(count,j) = nanmean(M_x(A,j));
    end
    clear A
end
meses = {'Ene','Feb','Mar','Abr','May','Jun','Jul','Ago','Sep','Oct','Nov','Dic'};

figure()
plot(dat(:,1),'o-','LineWidth',2)
hold on
plot(dat(:,2),'o-','LineWidth',2)
plot(dat(:,3),'o-','LineWidth',2)
xlim([1 12])
xticklabels(meses)
grid minor
xlabel('Meses')
ylabel('Transporte de Ekman promedio [m^3/s/m]')
title('Ciclo anual del transporte de Ekman')
legend('37°S','30°S','21°S',Location='northeast')

%%
ci = 2.5; % m/s

for i = 1:3
        Lr(i) = ci ./ f(i); % m
    w(:,i) = M_x(:,i)/Lr(i); % m/s
end

w = w * 86400; % m/dia
estaciones = {'Verano','Otoño','Invierno','Primavera'};

A =  [0 0 1; 1 0 0; 0.9290 0.6940 0.1250];

B = [2 3 4 1];
for k = 1:3
    figure()
    for i = 0:3
        A1 = find(fecha(:,2) == (3*i + 1));%
        A2 = find(fecha(:,2) == (3*i + 2));
        A3 = find(fecha(:,2) == (3*i + 3));
        est = [A1' A2' A3'];
        subplot (2,2,B(i+1))
            histogram(w(est,k),"NumBins", 15, "BinEdges", [-5:10],'FaceColor',A(k,:), 'Normalization', 'probability')
            ylabel('[%]','FontSize',20)
            xlabel('w [m/día]','FontSize',20) 
            title(estaciones(i+1),'FontSize',20)
            grid minor
            sgtitle("Velocidades verticales en " + lat(k) + "°S",'FontSize',25,'FontWeight','bold')
        
        clear A1 A2 A3 est
    end
end


%%

A1 = find(fecha(:,1) == 2000 & fecha(:,2) == 1);
A2 = find(fecha(:,1) == 2000 & fecha(:,2) == 2);
A3 = find(fecha(:,1) == 2000 & fecha(:,2) == 3);

ver = [A1',A2',A3'];

B = fecha(ver,:);
tau = tau_y(ver,1);

media = nanmean(tau);
clear A1 A2 A3

figure()
bar(datetime(B),tau)
hold on
%yline(media,Color='r',LineWidth=2)
axis  tight
grid minor
xlabel('Fecha','FontSize',20)
ylabel('Esfuerzo del viento [Pa]','FontSize',20)
title('Esfuerzo del viento paralelo a la costa durante el verano del 2000 en 37°S','FontSize',25)

%%

D1= find(B(:,2)==1 & B(:,3)>6 & B(:,3)<12);
D2 = find(B(:,2)==1 & B(:,3)>12 & B(:,3)<28);

%Impulso tiene unidades de m/s, por lo que la integral debe resultar kg/ms
date1 = B(D1,3)*86400; % seg
date2 = B(D2,3)*86400; % seg

int1 = trapz(date1, tau(D1)); % kg/ms
int2 = trapz(date2, tau(D2)); % kg/ms

H = 50; % m

I1 = int1/(rho_agua*H); % m/s
I2 = int2/(rho_agua*H); % m/s

c = 2.5; % m/s
g = (c^2)/H; %m/s^2

L_R = c/f(1); % m posicion vinal "x"
dist = linspace(0,L_R,1000); % m
R = c/f(1); % m

A1 = I1* sqrt(H/g); % m

d = (I2/f(1)) - R; % m
A2 = H*(exp(d/R)); % m

count = 0;
for k = 1:length(dist)
    count = count + 1;
    v1(count) = A1*(sqrt(g/H))*(exp((-dist(k))/R)); % m/s
    v2(count) = A2*(sqrt(g/H))*(exp((-dist(k))/R)); % m/s
end

figure()
plot(-dist/1000,v1,'LineWidth',2)
hold on
plot(-dist/1000,v2,'LineWidth',2)
%plot(-d/1000,2.5,'ok',LineWidth=3)
axis tight
grid minor
legend('1^{er} Evento','2^{do} Evento',Location='northeast')
xlabel('Distancia desde la costa [km]')
ylabel('Velocidad de la parcela [m/s]')
title('Perfil de chorro costero para dos eventos en verano del 2000')

%%

surg37 = w(:,1);

B = find(surg37 < 1.5);
B = [0 B'];

count = 0;
for i = 2: length(B)
    D(i-1) = (B(i)-1)-B(i-1);
    if D(i-1) == 0
        continue;
    else
        count = count + 1;
        F(count) = B(i-1)+1; %para fecha incial
        L(count) = B(i)-1; %para decha final
    end
end
D(find(D == 0)) = [];

date_i = fecha(F',1:3); 
date_f = fecha(L',1:3);
% primer evento = fecha(1)
% 2do evento B(i-1)+1

figure()
subplot(1,4,[1 2 3])
    plot(datetime(date_i),D,'+')
    xlim([datetime([1999 1 1]) datetime(fecha(end,:))])
    xlabel('Fecha de inicio del evento')
    ylabel('Cantidad de días')
    title('En función de fecha inicial')
    grid minor
subplot(1,4,4)
    histogram(D,"BinEdge",[0:25],'FaceColor','b','Normalization','probability')
    %axis tight
    grid minor
    xlim([0 25])
    xlabel('Cantidad de días')
    ylabel('Porcentaje [%]')
    title('Frecuencia de la duración')
sgtitle('Duración de eventos de surgencia en 37°S')

%%
for i = 1:388
T(i) = hours(datetime(date_f(i,:)) - datetime(date_i(i,:))+1)/24;
X1(i) = (find(fecha(:,1) == date_i(i,1) & fecha(:,2)==date_i(i,2) & fecha(:,3) ==date_i(i,3)));
X2(i) = find(fecha(:,1) == date_f(i,1) & fecha(:,2)==date_f(i,2) & fecha(:,3) ==date_f(i,3));
durac = [X1(i):X2(i)]*86400;

     if T(i) == 1
         integral(i) = tau_y(X1(i),1); % kg/ms
     else
         
         integral(i) = trapz(durac, tau_y(X1(i):X2(i),1)); % kg/ms
     end
% 
 H = 50; % m
 Impulso(i) = integral(i)/(rho_agua*H); % m/s
 clear durac
end
figure()
subplot(1,4,[1 2 3])
    plot(datetime(date_i),Impulso,'+')
    grid minor
    xlim([datetime([1999 1 1]) datetime(fecha(end,:))])
    xlabel('Fecha de inicio del evento')
    ylabel('Intensidad del evento [m/s]')
    title('En función del tiempo')
subplot(1,4,4)
    histogram(Impulso,"NumBins", 16,'FaceColor','b','Normalization','probability')
    xlabel('Intensidad del evento [m/s]')
    ylabel('Porcentaje [%]')
    title('Frecuencia de la Intensidad')
sgtitle('Intensidad de eventos de surgencia en 37°S')


