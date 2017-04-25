function [es_theta, es_phi] = mie_pec(radius, frequency, theta, phi, nmax)

% Compute the complex-value scattered electric far field of a perfectly
% conducting sphere using the mie series. Follows the treatment in
% Chapter 3 of 
%
% Ruck, et. al. "Radar Cross Section Handbook", Plenum Press, 1970.
%  
% The incident electric field is in the -z direction (theta = 0) and is
% theta-polarized. The time-harmonic convention exp(jwt) is assumed, and
% the Green's function is of the form exp(-jkr)/r.
% 
% Inputs:
%   radius: Radius of the sphere (meters)
%   frequency: Operating frequency (Hz)
%   theta: Scattered field theta angle (radians)
%   phi: Scattered field phi angle (radians)
%   nmax: Maximum mode number for computing Bessel functions
% Outputs:
%   es_theta: Theta-polarized electric field at the given scattering angles
%   es_phi: Phi-polarized electric field at the given scattering angles
%
%   Output electric field values are normalized such that the square of the
%   magnitude is the radar cross section (RCS) in square meters.
%
%   Author: Walton C. Gibson, Tripoint Industries, Inc.

% speed of light
c = 299792458.0;

% radian frequency
w = 2.0*pi*frequency;

% wavenumber
k = w/c;

% conversion factor between cartesian and spherical Bessel/Hankel function
s = sqrt(pi/2*k*radius); 

% mode numbers
mode = 1:nmax; 

% compute spherical bessel, hankel functions
[J1(mode)] = besselj(mode + 1/2, k*radius); J1 = J1*s;
[H1(mode)] = besselh(mode + 1/2, 2, k*radius); H1 = H1*s;
[J2(mode)] = besselj(mode + 1/2 - 1, k*radius); J2 = J2*s;
[H2(mode)] = besselh(mode + 1/2 - 1, 2, k*radius); H2 = H2*s;

% 
% if any(any(ierr))                
%     disp('Warning: There was an accuracy error in evaluating a Bessel or Hankel function.')
% end

% Ruck, et. al. (3.2-1)
A(mode) = -((-i).^mode) .* ( J1 ) ./ ( H1 ) .* (2*mode + 1) ./ (mode.*(mode + 1));

% Ruck, et. al. (3.2-2), using derivatives of bessel functions 
B(mode) = ((-i).^(mode+1)) .*(k*radius*J2 - mode .* J1 ) ./ (k*radius*H2 - mode .* H1 ).* (2*mode + 1) ./ (mode.*(mode + 1));
     
sintheta = sin(theta);
costheta = cos(theta);
 
% first two values of the Associated Legendre Polynomial
plm(1) = -sintheta;
plm(2) = -3.0*sintheta*costheta;

S1 = 0.0;
S2 = 0.0;

p = plm(1);

% compute coefficients for scattered electric far field
for i_mode = 1:nmax

        % derivative of associated Legendre Polynomial
    if abs(costheta) < 0.999999
        if i_mode == 1
            dp = costheta*plm(1)/sqrt(1.0 - costheta*costheta);
        else
            dp = (i_mode*costheta*plm(i_mode) - (i_mode + 1)*plm(i_mode - 1))/sqrt(1.0 - costheta*costheta);
        end
    end
     
    if abs(sintheta) > 1.0e-6
        term1 = A(i_mode)*p/sintheta;
        term2 = B(i_mode)*p/sintheta;
    end

    if costheta > 0.999999
        % Ruck, et. al. (3.1-12)
        val = ((-i)^(i_mode-1))*(i_mode*(i_mode+1)/2)*(A(i_mode) + i*B(i_mode));
        S1 = S1 + val;
        S2 = S2 + val;
    elseif costheta < -0.999999
        % Ruck, et. al. (3.1-14)
        val = ((i)^(i_mode-1))*(i_mode*(i_mode+1)/2)*(A(i_mode) - i*B(i_mode));
        S1 = S1 + val;
        S2 = S2 - val;
    else
        % Ruck, et. al. (3.1-6)
        S1 = S1 + ((-i)^(i_mode+1))*(term1 + i*B(i_mode)*dp);
        % Ruck, et. al. (3.1-7)
        S2 = S2 + ((-i)^(i_mode+1))*(A(i_mode)*dp + i*term2);
    end
    
    % recurrence relationship for next Associated Legendre Polynomial
    if i_mode > 1
        plm(i_mode + 1) = (2.0*i_mode + 1)*costheta*plm(i_mode)/i_mode - (i_mode + 1)*plm(i_mode - 1)/i_mode;
    end
    p = plm(i_mode + 1);
end

% complex-value scattered electric far field, Ruck, et. al. (3.1-5)
es_theta = S1*cos(phi);
es_phi = -S2*sin(phi);

% normalize electric field so square of magnitude is RCS in square meters
es_theta = es_theta*sqrt(4.0*pi)/k;
es_phi = es_phi*sqrt(4.0*pi)/k;

return