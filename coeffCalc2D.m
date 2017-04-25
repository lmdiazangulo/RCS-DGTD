function [ coeff ] = coeffCalc2D(N)

Np = (N+1)*(N+2)/2;
nId = nodeIndices(N,2);
alpha = zeros(Np,N+1,N+1,N+1);
coeff = []; % Coeff stores the [node scale x1pow x2pow x3pow]
for i=1:Np
    alpha(i,:,:,:) = alphaPol2D( Rpol(nId(i,1),N), ...
        Rpol(nId(i,2),N), Rpol(nId(i,3),N));
    for j=1:(N+1)
        for k=1:(N+1)
            for l=1:(N+1)
                if alpha(i,j,k,l) ~= 0
                    coeff = cat(1,coeff, [i,alpha(i,j,k,l),j-1,k-1,l-1]);
                end
            end
        end
    end
end

end