function [ RCS ] = ...
  bistaticRCS(fileName, sym_xy, sym_xz, intFreq, thP, phP)
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
    fprintf('Converting %s.rcs ===>\n', fileName);
    fprintf('===> %s.mat\n', fileName);
    RCS2MAT(fileName);
% end
fprintf('Loading file: %s.mat\n', fileName);
load(cat(2,fileName,'.mat'));
if mod(numel(time),2)==1
    tmp = time(1:(numel(time)-1));
    clear time;
    time = tmp;
    clear tmp;
end
cd ..
% -------------------------------------------------------------------------
constants;
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
% -------------------------------------------------------------------------
% Fourier transformation of terms. 
fprintf('FT processing ... ');
% Computes FT.
nC = Ns;
[ExIncFT] = getDTFT(time, ExInc(1:Ns), intFreq, nC);
[EyIncFT] = getDTFT(time, EyInc(1:Ns), intFreq, nC);
[EzIncFT] = getDTFT(time, EzInc(1:Ns), intFreq, nC);
% Stores field data for the interest frequency.
ENear.x = zeros(Nm,1); ENear.y = zeros(Nm,1); ENear.z = zeros(Nm,1);
HNear.x = zeros(Nm,1); HNear.y = zeros(Nm,1); HNear.z = zeros(Nm,1);
for k=1:Nm
    ENear.x(k) = getDTFT(time, Ex(k,1:Ns), intFreq, nC);
    ENear.y(k) = getDTFT(time, Ey(k,1:Ns), intFreq, nC);
    ENear.z(k) = getDTFT(time, Ez(k,1:Ns), intFreq, nC);
    HNear.x(k) = getDTFT(time, Hx(k,1:Ns), intFreq, nC);
    HNear.y(k) = getDTFT(time, Hy(k,1:Ns), intFreq, nC);
    HNear.z(k) = getDTFT(time, Hz(k,1:Ns), intFreq, nC); 
end
fprintf('OK\n');
% Computes Near to Far Field transformation.
fprintf('Computing Near to Far field transformation ... ');
[ ERad ] = near2FarField(N, Nb, K, Nm, ... 
            v, intFreq, ENear, HNear, sym_xy, sym_xz, thP, phP);
fprintf('OK\n');
% Computes RCS from ERad.
fprintf('Computing Radar Cross Section ...');
beta = 2 * pi * intFreq / c0;
ERadSq = ERad .* ERad;
EincSq = abs(ExIncFT).^2 + abs(EyIncFT).^2 + abs(EzIncFT).^2 ;
RCS = (ERadSq .* beta.^2) / (4*pi .* EincSq);
fprintf('OK\n');
% 
% figure();
% plotRCSMesh(N, p.x, p.y, p.z);
