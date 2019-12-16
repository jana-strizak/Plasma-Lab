%% Paschen Curve 
clc
clear all
close all


P = [0.46 0.57 0.81 6.8 9.5 12]; %data collected during lab 
V = [646 506 404 570 600 630];
%V = (b*P)/(c+log(P))
%ln(V) = ln(b*P/(c+ln(P)))
%ln(V) = ln(b*P) - ln(c+ln(P)) = ln(b) + ln(P) - ln(c+ln(P)) 
%set x = ln(P), y = ln(V), ln(b) = d
% y = d + x - ln(c + x)

%another way
% y-d = x-ln(c+x)
x = log(P);
y = log(V);
[fit,gof,output]=paschen(x,y); %this will also plot the data and fit
error = confint(fit);
fig = gcf;
fig.Color = [1 1 1];

%assigning constants for regular equation V = (b*P)/(c+log(P))
c =fit.c;
cerr = error(2,1)-fit.c;
derr = error(2,2)-fit.d;
b = exp(fit.d);
berr = sqrt(exp(fit.d)*derr);

%min voltage and pressure V = (b*P)/(c+log(P))
syms V(P) 
V(P) = b*P/(c+log(P));
V_p = diff(V,P); 
P_min = eval(solve(V_p == 0, P));
P_min_err = exp(gof.rmse);
V_min = eval(V(P_min));
V_min_err = exp(gof.rmse);

figure();%plotting the function 
p1 = fplot(V); hold on;
plot(P_min,V_min,'r*'); hold on;
line([0,P_min],[V_min,V_min],'Color','r','LineStyle','--');hold on;
line([P_min,P_min],[0,V_min],'Color','r','LineStyle','--');hold off;
xlim([0,8]); ylim([0 2000]);
xlabel('Pressure (Torr)'); ylabel('Voltage (V)');
title('Paschen Curve for Air');
legend('Fit','Min Point');
grid on;
fig = gcf;
fig.Color = [1 1 1];

%% Display  Paschen fits 
cprintf('*blue',   '          Paschen Fit Results \n');
table([b;c],[berr;cerr],'VariableNames',{'Value','Uncertainty'},'RowNames',{'B','C'})



%% Langmuir Probe 
close all
clear all
clc
load('Langmuir.mat');
figure();

%0.8Torr
Scale_ch1 = 5;
Shift_ch1 = -16.6;
Scale_ch2 = 0.5;
Shift_ch2 = 0.42;
I1 = -(Langmuir1(:,1)/(10*1000));%+Shift_ch1)*Scale_ch1;
V1 = (Langmuir1(:,2));%+Shift_ch2)*Scale_ch2;

[fit1 gof1,h1] = LangmuirFit(V1,I1,0,3,0.7431,0.3922,[0, 0, 1],[0, 0.4470, 0.7410]);
title('Langmuir Fit at 0.8 Torr');
error1 = confint(fit1);

%0.4Torr
Scale_ch1 = 5;
Shift_ch1 = -16.6;
Scale_ch2 = 0.5;
Shift_ch2 = 1.1;
I2 = -(Langmuir2(:,1)/(10*1000));%+Shift_ch1)*Scale_ch1;
V2 = (Langmuir2(:,2));%+Shift_ch2)*Scale_ch2;

[fit2 gof2,h2] = LangmuirFit(V2,I2,0,2,1,3,[0, 0.5, 0],[0.4660, 0.6740, 0.1880]);
title('Langmuir Fit at 0.4 Torr');
error2 = confint(fit2);


%0.25Torr
Scale_ch1 = 5;
Shift_ch1 = -16.6;
Scale_ch2 = 0.5;
Shift_ch2 = 1.1;
I3 = -(Langmuir3(:,1)/(10*1000));%+Shift_ch1)*Scale_ch1;
V3 = (Langmuir3(:,2));%+Shift_ch2)*Scale_ch2;

[fit3 gof3,h3] = LangmuirFit(V3,I3,0.8626,0.7150,0.8786,0.2575,[1, 0, 0],[0.8500, 0.3250, 0.0980]);
title('Langmuir Fit at 0.25 Torr');
error3 = confint(fit3);

fig = gcf;
fig.Color = [1 1 1];
legend([h1(2),h2(2),h3(2)],{Pressure1,Pressure2,Pressure3});
%% Display Langmuir fits 
cprintf('*blue',   '                                         Langmuir Fit Results \n');

a = [fit1.a,fit2.a,fit3.a];
aerr=[error1(2,1)-fit1.a,error2(2,1)-fit2.a,error3(2,1)-fit3.a];

b = [fit1.b,fit2.b,fit3.b];
berr =[error1(2,2)-fit1.b,error2(2,2)-fit2.b,error3(2,2)-fit3.b];

c = [fit1.c,fit2.c,fit3.c];
cerr =[error1(2,3)-fit1.c,error2(2,3)-fit2.c,error3(2,3)-fit3.c];

d = [fit1.d,fit2.d,fit3.d];
derr =[error1(2,4)-fit1.d,error2(2,4)-fit2.d,error3(2,4)-fit3.d];

table(a',aerr',b',berr',c',cerr',d',derr','VariableNames',{'a','aError','b','bError','c','cError','d','dError'},'RowNames',{Pressure1,Pressure2,Pressure3})

%% Derivation of plasma perameters 
close all
clc

%plasma temperature 
% b= -e/Te ----> Te = -e/b
qe = 1.60217662e-19;%in c
Kb = 8.617333262145e-5; %in ev
Te = 1./(b); %in eV
Te_err = (berr./b).*Te;

%Plasma Potential 
%the shift in the x is the plasma potential 
% c = Vp
Vp=c;
Vp_err = cerr;

 
%Plasma Density 
% Related to slope of exponential 
% a = e*n_e*SA*sqrt(Te/2pi*m_e)
area = 1; %area of probe
area_err = 0.5;

m_e = 9.10938356e-31;
syms n_e(SA,a_sym,temp)
n_e(SA,a_sym,temp) = (a_sym)/(qe*SA*sqrt(temp/(2*pi*m_e)));

density = [];
density_err = [];

density(1) = eval(n_e(area,a(1),Te(1)));
density_err(1) = PropError(n_e,[SA,a_sym,temp],[area,a(1),Te(1)],[area_err,aerr(1),Te_err(1)]);

density(2) = eval(n_e(1,a(2),Te(2)));
density_err(2) = PropError(n_e,[SA,a_sym,temp],[area,a(2),Te(2)],[area_err,aerr(2),Te_err(2)]);

density(3) = eval(n_e(1,a(3),Te(3)));
density_err(3) = PropError(n_e,[SA,a_sym,temp],[area,a(3),Te(3)],[area_err,aerr(3),Te_err(3)]);

%% Display Results 
cprintf('*blue',   '                                         Langmuir Parameters \n');
table(Te',Te_err',Vp',Vp_err',density',density_err','VariableNames',{'Te','TeError','Vp','VpError','Density','DensityError'},'RowNames',{Pressure1,Pressure2,Pressure3})