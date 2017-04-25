function [ RCS ] = RCSint(N, fq, a, Ar, J, M, EInc)
% Computes the scattered field for numel(fq) frequencies.
% Salvador Garcia's thesis p.139
% Output:
% RCS = Scattered field for all frequencies introduced with fq.
% Inputs:
% N = Order of function basis.
% fq = frequencies to be computed.
% a = geometric factor including r.
% Ar = surface area of elements.
% J and M = Electric and magnetic currents on nodes in frequency domain.
% NOTES:
% - rMod will be considered 1.
% - Works only for N == 1.

WRONG

if N > 1
    disp('WARNING, finalInt: Geometric integral is not ready for that order N');
end

numFreq = numel(fq); % Number of frequencies to be computed.
K = numel(Ar); % Number of elements.

Np = (N+1)*(N+2)/2; % Number of nodes.

c0 = 3e8; % Speed of light
mu0 = pi*4e-7; % Permeability
c0mu0 = c0*mu0;

RCS = zeros(1, numFreq);

for j=1:numFreq

    % Computes common factor of the integral
    % NOTE: a 4pi term could be simplified, left here for clarity.
    beta = 2*pi*fq(j)/c0;   
    cFactor = beta^2 / (4*pi* abs(EInc(j))^2 );

    RCS(j) = 0;
    for e=0:(K-1)
        
        % Computes geometric integrals for that element.
        g1 = geomInt(beta, [a(e+1,1) a(e+1,2) a(e+1,3)]); % First node
        g2 = geomInt(beta, [a(e+1,2) a(e+1,3) a(e+1,1)]); % Second node
        g3 = geomInt(beta, [a(e+1,3) a(e+1,1) a(e+1,2)]); % Third node

        % Computes integral.
        integral = g1*(c0mu0*J(Np*e+1,j)-M(Np*e+1,j)) + ...
                   g2*(c0mu0*J(Np*e+2,j)-M(Np*e+2,j)) + ...
                   g3*(c0mu0*J(Np*e+3,j)-M(Np*e+3,j));

        % Adds contribution of element e.
        RCS(j) = RCS(j) + 2*Ar(e+1)*integral;
        
    end
    RCS(j) = cFactor * abs(RCS(j))^2;
     
end