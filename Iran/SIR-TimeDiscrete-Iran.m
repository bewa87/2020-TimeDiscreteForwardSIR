clf; clc; close all; clear all;
pkg load statistics;

%
% Step 1: Process Data
%

% Step 1.1: Load Population Data

population_iran = 82360000; % Statista.com (2018)
N               = population_iran;

% Step 1.2: Load Infected Cases

[infected1,infected2] = textread('Iran-Confirmed.csv','%s;%f','delimiter',';');
infected2             = infected2(end:-1:1);
dates                 = (1:1:length(infected2))';

% Step 1.3: Load Recovered Cases

[recovered1,recovered2] = textread('Iran-Recovered.csv','%s;%f','delimiter',';');
recovered2              = recovered2(end:-1:1);

% Step 1.4: Load Dead Cases

[dead1, dead2] = textread('Iran-Dead.csv','%s;%f','delimiter',';');
dead2          = dead2(end:-1:1);

% Step 1.5: Process Data according to Preprint

R              = recovered2 + dead2;
I              = infected2 - R;
S              = N - I - R;

% Step 1.6: Plots of unprocessed and processed Data

figure(1)
hold off
plot(dates, infected2, 'color', 'black', '-@');
lt01 = title("Iran: People vs. Dates (Unprocessed)", "fontsize", 14);
xt01 = xlabel("t", "fontsize", 14);
yt01 = ylabel("People", "fontsize", 14);
hold on
plot(dates, R, 'color', 'blue', '-@');
legend({"Infected", "Recovered"}, 'location', 'northwest');

figure(2)
hold off
plot(dates, I, 'color', 'black', '-@');
lt02 = title("Iran: People vs. Dates (Processed)", "fontsize", 14);
xt02 = xlabel("t", "fontsize", 14);
yt02 = ylabel("People", "fontsize", 14);
hold on
plot(dates, R, 'color', 'blue', '-@');
legend({"Infected", "Recovered"}, 'location', 'northwest');

%
% Step 2: Parameter Identification
%

% Step 2.1: Parameter Idenfication for Alpha(t) and Beta

A_help = S(2:1:end) - S(1:1:(end-1));
B_help = I(2:1:end) - I(1:1:(end-1));
C_help = R(2:1:end) - R(1:1:(end-1));

x_param_L1 = fminsearch (@(x) ( sum(abs(A_help+x(1)*exp(-x(2)*dates(2:1:end)).*(I(2:1:end).*S(2:1:end))/N) ...
                           + abs(B_help-x(1)*exp(-x(2)*dates(2:1:end)).*(I(2:1:end).*S(2:1:end))/N+x(3)*I(2:1:end)) ...
                           + abs(C_help-x(3)*I(2:1:end))) ), [0.4;0.1;0.25]);

p = 2.90;
                           
x_param_L2 = fminsearch (@(x) ( sum( (abs(A_help+x(1)*exp(-x(2)*dates(2:1:end)).*(I(2:1:end).*S(2:1:end))/N)).^p ...
                           + (abs(B_help-x(1)*exp(-x(2)*dates(2:1:end)).*(I(2:1:end).*S(2:1:end))/N+x(3)*I(2:1:end))).^p ...
                           + (abs(C_help-x(3)*I(2:1:end))).^p ) ), [0.4;0.1;0.25]);

% Plot 2.2 Plots of Alpha(t) and Beta

alpha_data = N*(S(1:1:(end-1)) - S(2:1:end))./(I(2:1:end).*S(2:1:end));
beta_data  = (R(2:1:end)-R(1:1:(end-1)))./(I(2:1:end));

figure(3)
hold off
plot(dates(2:1:end),alpha_data,'color','black','-@');
lt03 = title("Iran: Alpha vs. Dates", "fontsize", 14);
xt03 = xlabel("t", "fontsize", 14);
yt03 = ylabel("Alpha", "fontsize", 14);
hold on
plot(dates(2:1:end), x_param_L2(1)*exp(-x_param_L2(2)*dates(1:1:(end-1))), 'color', 'red', '-@');
legend({"Alpha (Data)", "Alpha (Model)"}, 'location', 'northwest');
       
figure(4)
hold off
plot(dates(2:1:end),beta_data,'color','black','-@');
lt03 = title("Iran: Beta vs. Dates", "fontsize", 14);
xt03 = xlabel("t", "fontsize", 14);
yt03 = ylabel("Beta", "fontsize", 14);
hold on
plot(dates(2:1:end), x_param_L2(3)*ones(length(dates(2:1:end)),1), 'color', 'red', '-@');
legend({"Beta (Data)", "Beta (Model)"}, 'location', 'northwest');

%
% Step 3: Solve implicit time-discrete solution scheme
%
                           
% Step 3.1: Inputs

alpha01   = x_param_L2(1);
alpha02   = x_param_L2(2);
beta_mean = x_param_L2(3);
alpha     = alpha01*exp(-alpha02*(dates-1));

% Step 3.2: Time steps

Delta     = 1.0;
S_comp    = zeros(dates,1);
I_comp    = zeros(dates,1);
R_comp    = zeros(dates,1);
S_comp(1) = S(1);
I_comp(1) = I(1);
R_comp(1) = R(1);

% Step 3.3: Solve the implicit time-discrete solution scheme

for j = 1:1:(length(dates)-1)
  A = (1 + beta_mean*Delta)*alpha(j+1)*Delta;
  B = ((1 + beta_mean*Delta)*N - alpha(j+1)*Delta*(S_comp(j) + I_comp(j)))/2;
  
  I_comp(j+1) = -B/A + sqrt(B*B+N*I_comp(j)*A)/A;
  S_comp(j+1) = S_comp(j)/(1+alpha(j+1)*Delta*I_comp(j+1)/N);
  R_comp(j+1) = R_comp(j)+beta_mean*Delta*I_comp(j+1);
endfor

% Step 3.4: Plots

figure(5)
hold off
plot(dates, I, 'color', 'black', '-@');
lt02 = title("Iran: Infected vs. Dates", "fontsize", 14);
xt02 = xlabel("t", "fontsize", 14);
yt02 = ylabel("Infected", "fontsize", 14);
hold on
plot(dates, I_comp, 'color', 'red', '-@');
legend({"Infected (Data)", "Infected (Model)"}, 'location', 'northwest');

figure(6)
hold off
plot(dates, R, 'color', 'black', '-@');
lt02 = title("Iran: Recovered vs. Dates", "fontsize", 14);
xt02 = xlabel("t", "fontsize", 14);
yt02 = ylabel("Recovered", "fontsize", 14);
hold on
plot(dates, R_comp, 'color', 'red', '-@');
legend({"Recovered (Data)", "Recovered (Model)"}, 'location', 'northwest');