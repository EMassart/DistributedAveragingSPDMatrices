
% Run the tests of the paper
% Estelle M. Massart, Julien M. Hendrickx, P.-A. Absil, Extending a
% two-variable mean to a multivariable mean, Proceedings of ESANN 2016, pp.
% 101-106, 2016

% Author: E. Massart

clc;
clear all;
data.number = 10;
data.size = 20;
data.n_test = 100;   
data.str = 'result_papierEsann1_end.mat';
compare_algorithms( data );


clc;
clear all;
data.number = 30;
data.size = 20;
data.n_test = 100;
data.str = 'result_papierEsann2_end.mat';
compare_algorithms( data );


%plot the results : first test
clear all;
load 'result_papierEsann1_end.mat';
figure;
x_prov = 1:length(distM1);
col0 = [0 0 0];
col1 = [238 99 99]./255;
col2 = [135 206 235]./255;
x = data.number*(x_prov-ones(size(x_prov)));
semilogy(x,distM1,'o-','Color',col0);
hold on;
semilogy(x,distM2,'o-','Color',col1);
semilogy(x,distM3,'o-','Color',col2);
semilogy(x,distMMin1,'x-','Color',col0);
semilogy(x,distMMin2,'x-','Color',col1);
semilogy(x,distMMin3,'x-','Color',col2);
semilogy(x,distMMax1,'s-','Color',col0);
semilogy(x,distMMax2,'s-','Color',col1);
semilogy(x,distMMax3,'s-','Color',col2);
xlabel('Number of two-variable means evaluations');
ylabel('Error E_{rel}');
legend('Cyclic','Cyclic\_Random','Cyclic\_Cheap');


%plot the error wit respect to the CPU time
figure;
col0 = [0 0 0];
col1 = [238 99 99]./255;
col2 = [135 206 235]./255;
semilogy(timeM1,distM1,'o-','Color',col0);
hold on;
semilogy(timeM2,distM2,'o-','Color',col1);
semilogy(timeM3,distM3,'o-','Color',col2);
semilogy(timeM1,distMMin1,'x-','Color',col0);
semilogy(timeM2,distMMin2,'x-','Color',col1);
semilogy(timeM3,distMMin3,'x-','Color',col2);
semilogy(timeM1,distMMax1,'s-','Color',col0);
semilogy(timeM2,distMMax2,'s-','Color',col1);
semilogy(timeM3,distMMax3,'s-','Color',col2);
xlabel('CPU time');
ylabel('Error E_{rel}');
legend('Cyclic','Cyclic\_Random','Cyclic\_Cheap');

% %plot the results : second test
clear all;
load 'result_papierEsann2_end.mat';
figure;
x_prov = 1:length(distM1);
col0 = [0 0 0];
col1 = [238 99 99]./255;
col2 = [135 206 235]./255;
x = data.number*(x_prov-ones(size(x_prov)));
semilogy(x,distM1,'o-','Color',col0);
hold on;
semilogy(x,distM2,'o-','Color',col1);
semilogy(x,distM3,'o-','Color',col2);
semilogy(x,distMMin1,'x-','Color',col0);
semilogy(x,distMMin2,'x-','Color',col1);
semilogy(x,distMMin3,'x-','Color',col2);
semilogy(x,distMMax1,'s-','Color',col0);
semilogy(x,distMMax2,'s-','Color',col1);
semilogy(x,distMMax3,'s-','Color',col2);
xlabel('Number of two-variable means evaluations');
ylabel('Error E_{rel}');
legend('Cyclic','Cyclic\_Random','Cyclic\_Cheap');



%plot the error wit respect to the CPU time
figure;
col0 = [0 0 0];
col1 = [238 99 99]./255;
col2 = [135 206 235]./255;
semilogy(timeM1,distM1,'o-','Color',col0);
hold on;
semilogy(timeM2,distM2,'o-','Color',col1);
semilogy(timeM3,distM3,'o-','Color',col2);
legend('Cyclic','Cyclic\_Random','Cyclic\_Cheap');
semilogy(timeM1,distMMin1,'x-','Color',col0);
semilogy(timeM2,distMMin2,'x-','Color',col1);
semilogy(timeM3,distMMin3,'x-','Color',col2);
semilogy(timeM1,distMMax1,'s-','Color',col0);
semilogy(timeM2,distMMax2,'s-','Color',col1);
semilogy(timeM3,distMMax3,'s-','Color',col2);
xlabel('CPU time');
ylabel('Error E_{rel}');
legend('Cyclic','Cyclic\_Random','Cyclic\_Cheap');


