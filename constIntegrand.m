function [cjik] = constIntegrand(sd,wd,N,Nb,v)

Nfp = (N+1)*(N+2)/2;
Nbp = (Nb+1)*(Nb+2)/2;
Ncp = length(wd);
K   = numel(v.x)/Nbp;

% Initializes basis functions to be evaluated in cubature points.
caiwj = zeros(Ncp,Nfp); 
% Calculates the coefficients of the basis functions.
coeff = coeffCalc2D(N); 
for j=1:Ncp
    for c=1:size(coeff,1)
        i = coeff(c,1);
        const = coeff(c,2);
        % Evaluates alpha functions in dunavant cubature points.
        caiwj(j,i) = caiwj(j,i) + ...
         const * sd(j,1)^coeff(c,3) * sd(j,2)^coeff(c,4) * sd(j,3)^coeff(c,5);      
    end
    % Multiplies dunavant weights by the alpha functions.
    caiwj(j,:)=wd(j).*caiwj(j,:);
end
% Evaluates derivative of alpha basis wrt simplex i in cubature points.
coeffb = coeffCalc2D(Nb);
cdabji = zeros(Ncp, Nbp, 3); 
for i=1:3
    % dc := Derivative of alpha basis with respect to simplex coordinate i
    dc = dcoeff(coeffb,i); 
    for j=1:Nbp
        edc = []; % edc := dc to be evaluated.
        for c=1:size(dc,1);
            if dc(c,1)==j
                edc = cat(1,edc,dc(c,:));
            end
        end
        % Evaluates edc in all dunavant cubature points.
        cdabji(:,j,i) = evalCoeff(edc,sd);
    end
end
% This will store the modulus of the surface vector n evaluated in 
% dunavant points.
nmod = zeros(Ncp,K); 
for k=1:K
    %Takes vertices belonging to the element.
    xb = v.x((k-1)*Nbp+(1:Nbp));
    yb = v.y((k-1)*Nbp+(1:Nbp)); 
    zb = v.z((k-1)*Nbp+(1:Nbp));
    % Builds derivatives of coordinates vector wrt simplex coordinates.
    % The derivative is evaluated in dunavant cubature points.
    dx=zeros(Ncp,3); dy=zeros(Ncp,3); dz=zeros(Ncp,3);
    for i=1:3 
        for j=1:Nbp
            dx(:,i)=dx(:,i)+xb(j).*cdabji(:,j,i);
            dy(:,i)=dy(:,i)+yb(j).*cdabji(:,j,i);
            dz(:,i)=dz(:,i)+zb(j).*cdabji(:,j,i);
        end
    end
    % Calculates tangent vectors.
    t1x=dx(:,1)-dx(:,3); t1y=dy(:,1)-dy(:,3); t1z=dz(:,1)-dz(:,3);
    t2x=dx(:,2)-dx(:,3); t2y=dy(:,2)-dy(:,3); t2z=dz(:,2)-dz(:,3);  
    % Computes normal vectors from tangent vectors.
    nx=t1y(:).*t2z(:)-t1z(:).*t2y(:);
    ny=t1z(:).*t2x(:)-t1x(:).*t2z(:);
    nz=t1x(:).*t2y(:)-t1y(:).*t2x(:);
    % Computes modulus. NOTE: Modulus multiplied by wd should be equal to
    % the surface area.
    nmod(:,k) = sqrt(nx(:).^2 + ny(:).^2 + nz(:).^2)/2;
end
% Computes final form of constant integrand part.
cjik = zeros(Ncp,Nfp,K);
for i=1:Nfp
    for k=1:K
        cjik(:,i,k) = caiwj(:,i).*nmod(:,k);
    end
end