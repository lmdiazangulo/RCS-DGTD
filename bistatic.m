clear all;
close all;

addpath ../gidutils

addpath /home/luis/workspace/Cudg3d/models/tests-cylinder.gid
constants;

sim{1}.fileName = 'cyl-s2c5cm-o1-lf2';
sim{1}.color = 'b';
sim{2}.fileName = 'cyl-s2c5cm-o2-lf2';
sim{2}.color = 'k';


figure(1);
hold on;
for i = 1:numel(sim)
    numData = GiDGraphReader([sim{i}.fileName '.default.post.rcs.XZ.1.5e+09.grf']);
%     excData = GiDGraphReader([fileName '.post.exc.pw.grf']);
    plot(rad2deg(numData(:,1)), numData(:,2),sim{i}.color);
    grid on;
    xlabel('\theta');
    ylabel('dBm'); 

end 
% ylim([-50 10]);

% intFreq=1.5e9;
% sphereRadius=1;
% hold on;
% plotMie;
% Eth = 10.^(Eth_db/10);
% Ephi = 10.^(Ephi_db/10);
% E = sqrt(Eth.^2 + Ephi.^2);
% hold on;
% plot(phi, 10.*log10(E),'.-r');

% figure(2);
% plot(excData(:,1), excData(:,2),'.-b');
% [fq excFq] = getFFT(excData(:,1), excData(:,2)', 10000);
% figure(3);
% plot(fq, abs(excFq));