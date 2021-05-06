%%
% clear;
clc;
close all;
setPlotParameters;
%% test Crow instability
a1 = 0.0985;
G1 = 1.;
a2 = 0.0985;
G2 = -1.;
b = 1.;
Re = 2E30;
figure;
wavelen = [0.01:0.05:15]';
sigma = [];
for lam = wavelen'
    k = 2 * pi / lam;
    L = growthrate(a1, G1, a2, G2, b, Re, k, 0);
    Emax = max(real(eig(L, 'nobalance')));
    sigma = [sigma; Emax];
end
plot(wavelen, sigma, '-k');
hold on
axis([0 15 0 1])
title('Crow instability')
xlabel('\lambda')
ylabel('\sigma')
xticks([0:1:15])
%% validate the formula of Leweke 2016
a1 = 0.22;
G1 = 1.;
a2 = 0.22;
G2 = -1.;
b = 1.;
Re = 2000;
figure;
wavelen = [0.01:0.05:15]';
sigma = [];
for lam = wavelen'
    k = 2 * pi / lam;
    L = growthrate(a1, G1, a2, G2, b, Re, k, 1);
    Emax = max(real(eig(L, 'nobalance')));
    sigma = [sigma; Emax];
end
plot(wavelen, sigma, '-k');
hold on
axis([0 15 0 1])
title('Check formula')
xlabel('\lambda')
ylabel('\sigma')
%leweke
dataleweke = [4.160467587672689, 0.002034019913765839
4.176408076514347, 0.05139876943487409
4.192348565356004, 0.11987179921075974
4.240170031880977, 0.20586413694605943
4.3039319872476085, 0.2918581668776271
4.415515409139212, 0.3937808402769788
4.527098831030818, 0.4734105200457569
4.7502656748140275, 0.567383255379492
5.0371944739638685, 0.6518086193709092
5.260361317747078, 0.6964182973798032
5.595111583421891, 0.734670394010979
5.89798087141339, 0.7569955393706383
6.26461211477152, 0.7681809567000819
6.583421891604678, 0.77299187068913
6.918172157279491, 0.7698426934349554
7.157279489904356, 0.7666833630031746
7.380446333687566, 0.761929983687228
7.763018065887353, 0.750824099582366
8.257173219978746, 0.7365453474755816
8.75132837407014, 0.717489525305103
9.341126461211477, 0.6968514996243325
10.600425079702445, 0.6476221258046393
11.445270988310307, 0.6158646784488653
12.274176408076515, 0.5856978955847215
13.055260361317748, 0.5634878195712651
13.565356004250797, 0.5460260462849523
14.075451647183849, 0.5301566296865375
14.984059511158344, 0.5063677345553248];
plot(dataleweke(:,1), dataleweke(:,2))
%% Donnadieu 2009
a = 0.2066 * 1.36;
G = 1.;
b = 1.;
Re = 2000;
k = 0.2/a1;
L = growthrate(a, G, a, -G, b, Re, k, 1);
Emax = max(real(eig(L, 'nobalance')));
sigma =  Emax * G / (2*pi*b*b)
sigma * 2
growth = [2.5 7.9506e+00 1.3199e+00 5.8652e-06
5 2.0692e+01 9.6440e-01 2.3887e-05
10 9.6838e+01 7.2782e-01 1.0169e-06
15 2.9712e+02 6.0417e-01 5.8985e-09
20 7.4821e+02 5.2662e-01 1.5418e-10
25 1.6940e+03 4.7332e-01 4.3966e-10
30 3.7198e+03 4.3616e-01 7.8889e-10
35 8.1937e+03 4.0976e-01 3.9283e-09
40 1.8089e+04 3.9005e-01 5.2328e-09
45 3.9774e+04 3.7458e-01 1.1305e-08];
figure
time = 0.1 * growth(:, 1) * 2 * pi * b * b / G;
plot(growth(:, 1), 2 * (sigma * time), '-k');
hold on
plot(growth(:, 1), log(growth(:, 2)), '-r');
%% theoretical line
a1 = 0.2;
a2 = 0.2;
G1 = 1.74;
G2 = -0.73;
b = 0.4;
Re = G1*1000;
%two separate cases
wavelen = [0.6 1.3]';
sigma = [];
for lam = wavelen'
    k = 2 * pi / lam;
    L = growthrate(a1, G1, a2, G2, b, Re, k, 1);
    Emax = max(real(eig(L, 'nobalance')));
    sigma = [sigma; Emax * G1 / (2*pi*b*b)];
end
sigma * 2
%continuous line
figure;
wavelen = [0.01:0.05:5]';
sigma = [];
for lam = wavelen'
    k = 2 * pi / lam;
    L = growthrate(a1, G1, a2, G2, b, Re, k, 1);
    Emax = max(real(eig(L, 'nobalance')));
    sigma = [sigma; Emax * G1 / (2*pi*b*b)];
end
plot(wavelen, sigma, '-k');
hold on
%CFD res
time15 = [4 9.1224e+01 3.0089e+00 7.6758e-06
5 1.0437e+02 3.0987e+00 1.5138e-05
6 1.0894e+02 3.1272e+00 1.0701e-05
7 1.0870e+02 3.1257e+00 7.0534e-06
8 1.0539e+02 3.1051e+00 6.5535e-06
10 9.4363e+01 3.0314e+00 8.1253e-06
12 8.2505e+01 2.9419e+00 5.7281e-05
14 7.2439e+01 2.8552e+00 9.1498e-05
16 6.4599e+01 2.7788e+00 5.1211e-06
18 5.8665e+01 2.7146e+00 1.2278e-05
20 5.4191e+01 2.6617e+00 1.9697e-09
30 4.3319e+01 2.5124e+00 6.7726e-05
40 3.9970e+01 2.4588e+00 6.1412e-06
50 3.8889e+01 2.4405e+00 3.8895e-06];
time5 = [4 1.6981e+04 1.9480e+00 4.1659e-05
6 6.3381e+04 2.2114e+00 1.9761e-06
8 1.0933e+05 2.3204e+00 9.0298e-05
9 1.1559e+05 2.3316e+00 2.4990e-05
10 1.1121e+05 2.3238e+00 1.1201e-05
11 1.0031e+05 2.3032e+00 5.0591e-06
12 8.6691e+04 2.2740e+00 2.3003e-06
14 6.0261e+04 2.2013e+00 2.8042e-06
16 4.0297e+04 2.1208e+00 1.9815e-06
18 2.6922e+04 2.0401e+00 1.2567e-06
20 1.8308e+04 1.9630e+00 4.7573e-07
22 1.2780e+04 1.8911e+00 7.8842e-07
24 9.1843e+03 1.8250e+00 8.9043e-07
26 6.7952e+03 1.7648e+00 3.1133e-06
28 5.1704e+03 1.7101e+00 5.4233e-06
30 4.0398e+03 1.6608e+00 9.9052e-06
40 1.6551e+03 1.4823e+00 3.8657e-07
50 1.0287e+03 1.3872e+00 8.5744e-05];
wavelen = [4:2:50]';
rate = (log(interp1(time5( :,1), time5( :,2), wavelen)) - ...
        log(interp1(time15(:,1), time15(:,2), wavelen)) ) ...
           / 3.5 / 2 ;
plot(wavelen / 10., rate, 'b-')
legend('theory', 'present CFD', 'Location', 'Best')
title('growth rate of 2:1 vortex pair')
xlabel('\lambda')
ylabel('\sigma')
axis([0 5 0 2])
%% theoretical lines
figure;
G2 = -1.;
b = 0.4;

a1 = 0.2;
a2 = 0.2;
symbols = {'-.k', '-k', '--k'}
for G1=1:1:3
    Re = G1 * 1000;
    wavelen = [0.01:0.05:5]';
    sigma = [];
    for lam = wavelen'
        k = 2 * pi / lam;
        L = growthrate(a1, G1, a2, G2, b, Re, k, 1);
        Emax = max(real(eig(L, 'nobalance')));
        sigma = [sigma; Emax * G1 / (2*pi*b*b)];
    end
    plot(wavelen, sigma, symbols{G1});
    hold on
end

a1 = 0.16;
a2 = 0.16;
symbols = {'-.r', '-r', '--r'}
for G1=1:1:3
    Re = G1 * 1000;
    wavelen = [0.01:0.05:5]';
    sigma = [];
    for lam = wavelen'
        k = 2 * pi / lam;
        L = growthrate(a1, G1, a2, G2, b, Re, k, 1);
        Emax = max(real(eig(L, 'nobalance')));
        sigma = [sigma; Emax * G1 / (2*pi*b*b)];
    end
    plot(wavelen, sigma, symbols{G1});
    hold on
end
legend('\Gamma_{TEV}:\Gamma_{LEV} = 1:-1', '2:-1', '3:-1', 'Location', 'Best')
title('growth rate of vortex pair')
xlabel('\lambda')
ylabel('\sigma')
axis([0 5 0 3])