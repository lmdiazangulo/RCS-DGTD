function plotRCSMesh(N,x,y,z)

Np = (N+1)*(N+2)/2;
K = numel(x)/Np;

EToV=[];
for e=0:(K-1)
    for j=1:(N^2)
        EToV = cat(1,EToV, e*Np+triTesselator(j));
    end
end

trimesh( EToV, x,y,z );
colormap([0 0 0]);
axis equal;
drawnow;

end