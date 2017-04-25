clear all;
% close all;
% -------------------------------------------------------------------------
fprintf('-------------------------------------------------------------\n');
fprintf('----------  Radar Cross Section Analysis --------------------\n');
fprintf('----- L. D. Angulo @ University of Granada On May, 2012 -----\n');
fprintf('-------------------------------------------------------------\n');
% -------------------------------------------------------------------------
addpath Dunavant
% -------------------------------------------------------------------------
constants;
% -------------------------------------------------------------------------
fileName = 'rcs-q20cm-o2-rk.default.post';
color = '.-b';
sym_xy = 2; % 0 for no symmetry, 1 for PEC, 2 for PMC.
sym_xz = 1; % 0 for no symmetry, 1 for PEC, 2 for PMC.
% -------------------------------------------------------------------------
% Converts .rcs to .mat
cd Data
% If fileName.mat does not exist it is converted to mat format. 
% exists = (exist([fileName '.mat'],'file')~=0);
% % Checks if file is outdated.
% if (exists)
%     RCSfile = dir([fileName '.rcs']);
%     RCSfiledate = datenum(RCSfile.date);
%     MATfile = dir([fileName '.mat']);
%     MATfiledate = datenum(MATfile.date);
%     outdated = (RCSfiledate > MATfiledate);
% end
% % If the .mat file does not exists or is outdated then rcs is converted.
% if (~exists || outdated)
%     fprintf('%s.mat is outdated. \n', fileName);
%     fprintf('It will be converted again.\n');
%     fprintf('Converting %s.rcs ===>\n', fileName);
%     fprintf('===> %s.mat\n', fileName);
    RCS2MAT(fileName);
% end
fprintf('Loading file: %s.mat\n', fileName);
load(cat(2,fileName,'.mat'));
% if mod(numel(time),2)==1
%     tmp = time(1:(numel(time)-1));
%     clear time;
%     time = tmp;
%     clear tmp;
% end
cd ..
% -------------------------------------------------------------------------
% Field data is assumed to be ordered with nodes in rows. Different columns
% represent different values of those nodes.
% Initializes variables.
Nm = numel(x);       % Number of map points.
Ns = numel(time);    % Number of samples.
Np = (N+1)*(N+2)/2;  % Number of nodes.
Nfp = (N+1)*(N+2)/2; % Number of face points
K = Nm/Nfp;          % Number of face elements
Nbp = size(v.x,1)/K; % Number of element points.
% Order of elements.
switch Nbp 
    case 3,  Nb = 1; 
    case 6,  Nb = 2; 
end
p.x = x;
p.y = y;
p.z = z;
%% ------------------------------------------------------------------------
intFreq  = 150e6;         % Interest frequency
numzeros = 10000;           % num of zeros to add for the FFT
thP = pi/2;  % Stores theta components.
phP = linspace(0, pi, 200);  % Stores phi components.
sphereRadius = 1;
% Fourier transformation of terms. 
% TODO This can be improved using DTFT (for bistatic RCS) rather than FFT.
% Extracting sampling time as an average of all dt.
sTime = time(2:numel(time)) - time(1:(numel(time)-1));
sTime = sum(sTime)/numel(sTime);
fprintf('Maximum frequency for this data is: %.3e [Hz]\n', 1/(sTime*2));
% Computes FFT.
% [fq ExIncFT] = getFFT(sTime, ExInc(1:Ns),numzeros);
% [fq EyIncFT] = getFFT(sTime, EyInc(1:Ns),numzeros);
% [fq EzIncFT] = getFFT(sTime, EzInc(1:Ns),numzeros);

% Obtains nearest frequency to the interest frequency.
% indexInt := index of interest.
% indexInt = find(min(abs(fq-intFreq))==abs(fq-intFreq));
% availIntFreq = fq(indexInt);
% fprintf('Interest frequency was: %.5e [Hz]\n', intFreq);
% fprintf('Nearest frequency found after FFT was: %.5e [Hz]\n',availIntFreq);
% Stores incident field data.
% ExIncIntFq = ExIncFT(indexInt);
% EyIncIntFq = EyIncFT(indexInt);
% EzIncIntFq = EzIncFT(indexInt);
% Stores field data for the interest frequency.
ENear.x = zeros(Nm,1); ENear.y = zeros(Nm,1); ENear.z = zeros(Nm,1);
HNear.x = zeros(Nm,1); HNear.y = zeros(Nm,1); HNear.z = zeros(Nm,1);
DTFTENear.x = zeros(Nm,1); DTFTENear.y = zeros(Nm,1); DTFTENear.z = zeros(Nm,1);
DTFTHNear.x = zeros(Nm,1); DTFTHNear.y = zeros(Nm,1); DTFTHNear.z = zeros(Nm,1);
for k=1:Nm
%     [fq ENear.x(k)] = getFFTintFq(sTime, Ex(k,:),numzeros, indexInt);
%     [fq ENear.y(k)] = getFFTintFq(sTime, Ey(k,:),numzeros, indexInt);
%     [fq ENear.z(k)] = getFFTintFq(sTime, Ez(k,:),numzeros, indexInt);
%     [fq HNear.x(k)] = getFFTintFq(sTime, Hx(k,:),numzeros, indexInt);
%     [fq HNear.y(k)] = getFFTintFq(sTime, Hy(k,:),numzeros, indexInt);
%     [fq HNear.z(k)] = getFFTintFq(sTime, Hz(k,:),numzeros, indexInt); 
    [DTFTENear.x(k)] = getDTFT(intFreq, time, Ex(k,:));
    [DTFTENear.y(k)] = getDTFT(intFreq, time, Ey(k,:));
    [DTFTENear.z(k)] = getDTFT(intFreq, time, Ez(k,:));
    [DTFTHNear.x(k)] = getDTFT(intFreq, time, Hx(k,:));
    [DTFTHNear.y(k)] = getDTFT(intFreq, time, Hy(k,:));
    [DTFTHNear.z(k)] = getDTFT(intFreq, time, Hz(k,:));
end
[ExIncDTFT] = getDTFT(intFreq, time, ExInc);
[EyIncDTFT] = getDTFT(intFreq, time, EyInc);
[EzIncDTFT] = getDTFT(intFreq, time, EzInc);
fprintf('FFT processing has finished.\n');
% Computes Near to Far Field transformation.
fprintf('Computing Near to Far field transformation... ');
% [ ERad ] = near2FarField(N, Nb, K, Nm, ...
%             v, availIntFreq, ENear, HNear, sym_xy, sym_xz, thP, phP);
[ ERad ] = near2FarField(N, Nb, K, Nm, ...
    v, intFreq, DTFTENear, DTFTHNear, sym_xy, sym_xz, thP, phP);
fprintf('OK\n');
% Computes RCS from ERad.
fprintf('Computing Radar Cross Section ...');
% beta = 2 * pi * availIntFreq / c0;
ERadSq = ERad .* ERad;
% EincSq = abs(ExIncIntFq).^2 + abs(EyIncIntFq).^2 + abs(EzIncIntFq).^2 ;
beta = 2 * pi * intFreq / c0; 
EincSq = abs(ExIncDTFT).^2 + abs(EyIncDTFT).^2 + abs(EzIncDTFT).^2 ;
RCS = (ERadSq .* beta.^2) / (4*pi .* EincSq);
fprintf('OK\n');
% -------------------------------------------------------------------------

fprintf('Starting plotting routines\n');
figure(1);
hold on;
plot(phP, 10*log10(RCS(1,:)),color);
xlim([0 pi]); 
xlabel('\phi');
ylabel('dBm'); 
grid on;
title([fileName ' @ ' num2str(intFreq*1e-6) ' [MHz]']);
% toc
fprintf('----------------------- END ---------------------------------\n');
plotMie;
