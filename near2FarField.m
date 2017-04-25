function [ ERad ] = ...
 near2FarField(N, Nb, K, Nm, v, fq, ENear, HNear, sym_xy, sym_xz, thP, phP)
% -------------------------------------------------------------------------
Nc = 10;             % Order of cubature integration rule
% -------------------------------------------------------------------------
% Some constants.
constants;
% -------------------------------------------------------------------------
Nfp = (N+1)*(N+2)/2; % Number of face points
Nbp = (Nb+1)*(Nb+2)/2; % Number of face points of the geometrical base.
% Calculates cubature parameters.
[xyd, wd]=dunavant_rule(Nc);
% Converts to simplex coordinates.
sd = [1-xyd(1,:)'-xyd(2,:)', xyd(1,:)', xyd(2,:)']; 
% Number of cubature points.
Ncp = length(wd);
% Computes J and M currents in surfaces. It is considered the fact that
% possibly n is not constant through the surface of the element.

% Does .rcs file contains symmetries?
% if sym_xy == 1 % PEC in plane XY
%     v.x = [ v.x;  v.x];
%     v.y = [ v.y;  v.y]; 
%     v.z = [ v.z; -v.z];
%     HNear.x = [ HNear.x, -HNear.x];
%     HNear.y = [ HNear.y, -HNear.y];
%     HNear.z = [ HNear.z,  HNear.z];
%     ENear.x = [ ENear.x,  ENear.x];
%     ENear.y = [ ENear.y,  ENear.y];
%     ENear.z = [ ENear.z, -ENear.z];
%     K  = K  * 2;
%     Nm = Nm * 2;
% end
if sym_xy == 2 % PMC in plane XY
    v.x = [ v.x;  v.x];
    v.y = [ v.y;  v.y]; 
    v.z = [ v.z; -v.z];
    HNear.x = [ HNear.x,  HNear.x];
    HNear.y = [ HNear.y,  HNear.y];
    HNear.z = [ HNear.z, -HNear.z];
    ENear.x = [ ENear.x, -ENear.x];
    ENear.y = [ ENear.y, -ENear.y];
    ENear.z = [ ENear.z,  ENear.z];
    K  = K  * 2;
    Nm = Nm * 2;
end
if sym_xz == 1 % PEC in plane XZ
    v.x = [ v.x;  v.x];
    v.y = [ v.y; -v.y];
    v.z = [ v.z;  v.z];
    HNear.x = [ HNear.x, -HNear.x];
    HNear.y = [ HNear.y,  HNear.y];
    HNear.z = [ HNear.z, -HNear.z];
    ENear.x = [ ENear.x,  ENear.x];
    ENear.y = [ ENear.y, -ENear.y];
    ENear.z = [ ENear.z,  ENear.z];
    K  = K  * 2;
    Nm = Nm * 2;
end
% if sym_xz == 2 % PMC in plane XZ
%     v.x = [ v.x;  v.x];
%     v.y = [ v.y; -v.y];
%     v.z = [ v.z;  v.z];
%     HNear.x = [ HNear.x,  HNear.x];
%     HNear.y = [ HNear.y, -HNear.y];
%     HNear.z = [ HNear.z,  HNear.z];
%     ENear.x = [ ENear.x, -ENear.x];
%     ENear.y = [ ENear.y,  ENear.y];
%     ENear.z = [ ENear.z, -ENear.z];
%     K  = K  * 2;
%     Nm = Nm * 2;
% end

[J M] = buildCurrents(sd, wd, N, Nb, v, ENear, HNear);
% -------------------------------------------------------------------------
% cjik stores the integral part which is not in the exponential, evaluated
% in cubature points.
% cjik = zeros(Ncp,Nfp,K);
temp = constIntegrand(sd, wd, N, Nb, v);
cjik = zeros(Ncp,Nm);
for k=1:K
    for i=1:Nfp
        cjik(:,(k-1)*Nfp+i) = temp(:,i,k);
    end
end
% Initializes basis functions to be evaluated in cubature points.
cabij = zeros(Ncp,Nbp);
% Calculates the coefficients of the basis functions.
coeff = coeffCalc2D(Nb);
for j=1:Ncp
    for c=1:size(coeff,1)
        i = coeff(c,1);
        const = coeff(c,2);
        cabij(j,i) = cabij(j,i) + const ...
         * sd(j,1)^coeff(c,3) * sd(j,2)^coeff(c,4) * sd(j,3)^coeff(c,5);
    end
end
% - Starts calculations of integrals.
% Calculates integral term, needs to be calculated for all the
% frequencies of interest.
beta = 2 * pi * fq / c0;
% Builds geometric vector (a1, a2, a3). Changes cartesian coordinates to
% simplex coordinates and multiplies by (rx, ry, rz).
E.x = zeros(numel(thP), numel(phP));
E.y = zeros(numel(thP), numel(phP));
E.z = zeros(numel(thP), numel(phP));
for p=1:numel(thP)
    for q=1:numel(phP)
        % Position vector cartesians coordinates components
        pos.x = sin(thP(p)) * cos(phP(q));
        pos.y = sin(thP(p)) * sin(phP(q));
        pos.z = cos(thP(p));
        % Computes  cdijk, The exponential part without multiplying by
        % beta evaluated in all dunavant cubature points.
        cdjk=zeros(Ncp,K);
        for k=1:K
            for i=1:Nbp
                xb = v.x((k-1)*Nbp + i);
                yb = v.y((k-1)*Nbp + i); 
                zb = v.z((k-1)*Nbp + i); 
                cdjk(:,k) = cdjk(:,k) + ...
                 (pos.x * xb + pos.y * yb + pos.z * zb) .* cabij(:,i);
            end
        end       
        % Calculates nodes integrand values
        [int] = setIntegrand(pos, M, J);
        % Calculates cos and sin parts.
        cPart = zeros(Ncp,Nm); 
        sPart = zeros(Ncp,Nm);
        for k=1:K
            cTemp = cos(beta.*cdjk(:,k));
            sTemp = sin(beta.*cdjk(:,k));
            for i=1:Nfp
                cPart(:,(k-1)*Nfp+i) = cTemp;
                sPart(:,(k-1)*Nfp+i) = sTemp;
            end
        end
        phase = cPart + 1i * sPart;
        % Calculates Scatt field for all frequencies.
        temp = cjik.*phase;
        E.x(p,q) = sum(sum(temp.*int.x));
        E.y(p,q) = sum(sum(temp.*int.y));
        E.z(p,q) = sum(sum(temp.*int.z));
    end
end
% Calculates Erad.
ERadSq = abs(E.x).^2 + abs(E.y).^2 + abs(E.z).^2;
ERad   = sqrt(ERadSq);
