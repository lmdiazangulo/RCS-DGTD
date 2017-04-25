function [ res ] = triTesselator(tr)
% function [ v1 v2 v3 ] = triTesselator(tr)
% Divides a triangle with nodes 1,...,Np in triangles with vertices in all
% nodes.

N=1; Ntr=1;
while tr/Ntr>1
    N=N+1;
    Ntr=N*N;
end

fNode = N*(N-1)/2 + 1;
sNode = N*(N+1)/2 + 1;
if mod(N+tr,2)==0
    v1 = fNode + floor((tr-(N-1)^2)/2);
    v2 = sNode + floor((tr-(N-1)^2)/2);
    v3 = v2 + 1;
else
    v3 = sNode + floor((tr-(N-1)^2)/2);
    v2 = fNode + floor((tr-(N-1)^2)/2);
    v1 = v2-1;
end

res = [v1 v2 v3];
end