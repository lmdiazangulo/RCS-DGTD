function [J M] = buildCurrents(sd, wd, N, Nb, v, E, H)
% -------------------------------------------------------------------------
Nfp = (N+1)*(N+2)/2;
Nbp = (Nb+1)*(Nb+2)/2;
Ncp = length(wd);
K = size(v.x, 1)/Nbp;
Nm = K*Nfp;
% -------------------------------------------------------------------------
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
% -------------------------------------------------------------------------
J.x = zeros(Ncp,Nm); 
J.y = zeros(Ncp,Nm); 
J.z = zeros(Ncp,Nm);
M.x = zeros(Ncp,Nm); 
M.y = zeros(Ncp,Nm); 
M.z = zeros(Ncp,Nm);
for k=1:K
    % Takes vertices belonging to the element.
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
    t1x = dx(:,1) - dx(:,3);
    t1y = dy(:,1) - dy(:,3); 
    t1z = dz(:,1) - dz(:,3);
    t2x = dx(:,2) - dx(:,3);
    t2y = dy(:,2) - dy(:,3);
    t2z = dz(:,2) - dz(:,3);
    % Computes normal vectors from tangent vectors.
    nx = t1y.*t2z - t1z.*t2y;
    ny = t1z.*t2x - t1x.*t2z;
    nz = t1x.*t2y - t1y.*t2x;
    % Computes modulus. NOTE: Modulus multiplied by wd should be equal to
    % the surface area.
    nmod = sqrt(nx.^2 + ny.^2 + nz.^2 );
    % Normalizes vectors.
    nx = nx./nmod;
    ny = ny./nmod;
    nz = nz./nmod;
    % Builds currents in cubature points.
    for p = 1:Nfp
        % Points to component to be calculated.
        i = (k-1)*Nfp+p; 
        % Magnetic currents int, M = - n vecProd E
        M.x(:,i) = -ny(:) .* E.z(i) + nz(:) .* E.y(i);
        M.y(:,i) = -nz(:) .* E.x(i) + nx(:) .* E.z(i);
        M.z(:,i) = -nx(:) .* E.y(i) + ny(:) .* E.x(i);
        % Electric currents terms, J = n vecProd H
        J.x(:,i) =  ny(:) .* H.z(i) - nz(:) .* H.y(i);
        J.y(:,i) =  nz(:) .* H.x(i) - nx(:) .* H.z(i);
        J.z(:,i) =  nx(:) .* H.y(i) - ny(:) .* H.x(i);
    end 
end
