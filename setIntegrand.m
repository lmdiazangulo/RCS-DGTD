function [ int ] = setIntegrand(pos, M, J)

    c0 = 299792458.0; % Speed of light
    mu0 = pi*4e-7; % Permeability
    c0mu0 = c0*mu0;
   
    tmp.x = c0mu0 * (J.y * pos.z - J.z * pos.y) + M.x;
    tmp.y = c0mu0 * (J.z * pos.x - J.x * pos.z) + M.y;
    tmp.z = c0mu0 * (J.x * pos.y - J.y * pos.x) + M.z;
    
    int.x = tmp.y * pos.z - tmp.z * pos.y;
    int.y = tmp.z * pos.x - tmp.x * pos.z;
    int.z = tmp.x * pos.y - tmp.y * pos.x;
    
end