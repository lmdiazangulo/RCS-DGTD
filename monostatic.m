clear all;
close all;

addpath ../gidutils

addpath /home/luis/workspace/Cudg3d/models/tests-pml-rcs.gid

sim{1}.fileName = 'rcs-q20cm-o2-rk.default.post.rcs.monostatic.grf';
sim{1}.color = '-b';
sim{2}.fileName = 'pml-rcs-q20cm-o2-rk.default.post.rcs.monostatic.grf';
sim{2}.color = '-k';

figure(1);
xlabel('Frequency [Hz]');
ylabel('dBm'); 

nSim = numel(sim);
for i=1:nSim
    numData = GiDGraphReader(sim{i}.fileName);
    semilogx(numData(:,1), numData(:,2), sim{i}.color);
    hold on;
end



thP = 0;
phP = 0;

sphereRadius = 1;

nFreq = 400;
minFq = 10e6;
maxFq = 2000e6;
% freq = linspace(minFq, maxFq, nFreq);
freq = logspace(log10(minFq), log10(maxFq), nFreq);
RCSMie = zeros(nFreq, 1);
for i=1:nFreq
    scaledFreq = freq(i) * sphereRadius;
    [es_theta, es_phi] = mie_pec(1, scaledFreq, thP, phP, 50);
    RCSMie(i) = abs(es_theta^2 + es_phi^2) *  sphereRadius^2;
end

hold on;
semilogx(freq, 10*log10(RCSMie(:)),'.-r');
grid on;
hold off;

xlabel('Frequency [Hz]');
ylabel('RCS [dBm]');
xlim([minFq maxFq]);
ylim([-5 12]);