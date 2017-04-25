function RCSreader(fileName)

RCSGlobals;

% Storing raw data in fdata variable.
fid = fopen(fileName,'r'); % Opens file
fdata = fscanf(fid,'%c'); % Reads file and puts all in a string called fdata.
fclose(fid); % Closes file

% --- Extracts information.
% Extracts N
label = 'N:';
fpointer = findstr(label,fdata);
labelSize = numel(label);
fpointer = fpointer+labelSize+1;
label = 'NODETOL:';
lpointer = findstr(label,fdata)-1;
numdata = fdata(fpointer:lpointer);
N = str2num(numdata);

% Extracts NODETOL
label = 'NODETOL:';
fpointer = findstr(label,fdata);
labelSize = numel(label);
fpointer = fpointer+labelSize+1;
label = 'Nm';
lpointer = findstr(label,fdata)-1;
numdata = fdata(fpointer:lpointer);
NODETOL = str2num(numdata);

% Extracts nx ny nz
label = 'nx ny nz:';
fpointer = findstr(label, fdata);
labelSize = numel(label);
fpointer = fpointer+labelSize+1;
label = 'x y z:';
lpointer = findstr(label,fdata)-1;
numdata=fdata(fpointer:lpointer);
numdata=str2num(numdata);
Nv = size(numdata,1); %Number of vectors
nVect = reshape(numdata,Nv,3);
nx = nVect(1:Nv, 1);
ny = nVect(1:Nv, 2);
nz = nVect(1:Nv, 3);

    
% Extracts x y z
containsV=1;
label = 'x y z:';
fpointer = findstr(label, fdata);
labelSize = numel(label);
fpointer = fpointer+labelSize+1;
label = 'vx vy vz:';
lpointer = findstr(label,fdata)-1;
if numel(lpointer)==0
    containsV=0;
    label = 'END HEADER';
    lpointer = findstr(label,fdata)-1;
end
numdata = fdata(fpointer:lpointer);
numdata=str2num(numdata);
Nv = size(numdata,1); % Number of coordinates
nVect = reshape(numdata,Nv,3);
x = nVect(1:Nv,1);
y = nVect(1:Nv,2);
z = nVect(1:Nv,3);

% Extracts vx vy vz
if containsV==1
    label = 'vx vy vz:';
    fpointer = findstr(label, fdata);
    labelSize = numel(label);
    fpointer = fpointer+labelSize+1;
    label = 'END HEADER';
    lpointer = findstr(label,fdata)-1;
    numdata = fdata(fpointer:lpointer);
    numdata = str2num(numdata);
    Nv = size(numdata,1); % Number of coordinates
    nVect = reshape(numdata,Nv,3);
    v.x = nVect(1:Nv,1);
    v.y = nVect(1:Nv,2);
    v.z = nVect(1:Nv,3);
end

% Counts the number of steps saved in the file.
iStepPos = findstr('RCSSTEP:',fdata);
fStepPos = findstr('END RCSSTEP',fdata);
Nsteps = numel(iStepPos);
Nv = length(x);

time = zeros(1,Nsteps);
ExInc = zeros(1,Nsteps); EyInc = ExInc; EzInc = EyInc;
Ex = zeros(Nv,Nsteps); Ey = zeros(Nv,Nsteps); Ez = zeros(Nv,Nsteps);
Hx = zeros(Nv,Nsteps); Hy = zeros(Nv,Nsteps); Hz = zeros(Nv,Nsteps);

for i=1:Nsteps
    % Separates step data from the rest of data. and performs same
    % procedure as before but using stepData rather than fData.
    stepData = fdata(iStepPos(i):fStepPos(i));
    
    % Extracts time
    label = 'time:';
    fpointer = findstr(label, stepData);
    labelSize = numel(label);
    fpointer = fpointer+labelSize+1;
    label = 'ExInc EyInc EzInc:';
    lpointer = findstr(label,stepData)-1;
    numdata = stepData(fpointer:lpointer);
    time(i) = str2num(numdata);
    
    % Extracts ExInc, ...
    label = 'ExInc EyInc EzInc:';
    fpointer = findstr(label,stepData);
    labelSize = numel(label);
    fpointer = fpointer+labelSize+1;
    label = 'Ex Ey Ez:';
    lpointer = findstr(label,stepData)-1;
    numdata = stepData(fpointer:lpointer);
    Einc = str2num(numdata);
    ExInc(i)=Einc(1); EyInc(i)=Einc(2); EzInc(i)=Einc(3);
    
    % Extracts Ex,Ey,Ez
    label = 'Ex Ey Ez:';
    fpointer = findstr(label, stepData);
    labelSize = numel(label);
    fpointer = fpointer+labelSize+1;
    label = 'Hx Hy Hz:';
    lpointer = findstr(label,stepData)-1;
    numdata = stepData(fpointer:lpointer);
    numdata = str2num(numdata);
    Ex(:,i) = numdata(:,1);
    Ey(:,i) = numdata(:,2);
    Ez(:,i) = numdata(:,3);
    
    % Extracts Hx,Hy,Hz
    label = 'Hx Hy Hz:';
    fpointer = findstr(label, stepData);
    labelSize = numel(label);
    fpointer = fpointer+labelSize+1;
    lpointer = length(stepData)-1;
    numdata = stepData(fpointer:lpointer);
    numdata = str2num(numdata);
    Hx(:,i) = numdata(:,1);
    Hy(:,i) = numdata(:,2);
    Hz(:,i) = numdata(:,3);

    clear stepData;
    
end




