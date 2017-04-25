clear all;
close all;

addpath ../gidutils
%% Bistatic 1 GHz
addpath /home/luis/mirror/almond.gid/
fileName = 'lf2-o2.default.post.rcs.1e+09.grf';
numData = GiDGraphReader(fileName);

figure(1);
plot(rad2deg(numData(:,1)), numData(:,2),'.-b');
grid on;
xlabel('\phi');
ylabel('dBm'); 
ylim([-40 20]);

addpath ./nasaAlmondData
fileName = 'bies_pec_lfdg_1G.dat';
[phi Eth_db Eth_deg Ephi_db Ephi_deg] = nasaDataReader(fileName,2,182);
Eth = 10.^(Eth_db/10);
Ephi = 10.^(Ephi_db/10);
E = sqrt(Eth.^2 + Ephi.^2);
hold on;
plot(phi, 10.*log10(E),'.-r');


%% Monostatic 500 MHz - 2 GHz
addpath /home/luis/mirror/almond.gid/
fileName = 'lf2-o2.default.post.rcs.monostatic.grf';
numData = GiDGraphReader(fileName);

figure(2);
semilogx(numData(:,1), numData(:,2),'.-b');
grid on;
xlabel('Frequency [Hz]');
ylabel('dBm'); 
xlim([500e6 2e9])