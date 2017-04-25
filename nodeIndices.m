function [ nId ] = nodeIndices(N,D)
% Generates number of node Indices for order N and dimension D

switch D
    case 1
        
        Np = N+1;
        nId = zeros(Np,D+1);
        
        for i=0:N
            nId(i+1,:) = [N-i, i];
        end
    
    case 2
             
        Nfp = N+1;
        Np = (N+1)*(N+2)/2;
        
        P = Pmatrix2D(N);
        % numerates nodes
% % %         nId1 = [];
% % %         for i = 1:Nfp
% % %             temp(1:i) = N+1-i;
% % %             nId1 = cat(2, nId1, temp);
% % %         end
        for i = 1:Nfp
            nId1((1+i*(i-1)/2):Np) = N+1-i;
        end
        
        % P is the permutation matrix (same as in Dorig{l})
        nId = [nId1' P'*nId1' P*nId1'];

    case 3
        
        Np = (N+1)*(N+2)*(N+3)/6;
       % Nfp = (N+1)*(N+2)/2;
        
        % Numerates nodes
% % %         nId1 = [];
% % %         for i = 0:N
% % %             temp = ones(1, (i+1)*(i+2)/2 )*(N-i);
% % %             nId1 = cat(2, nId1, temp);
% % %         end
        for i=0:N
            nId1( (1+(i+2)*i*(i+1)/6):Np ) = N-i;
        end

        % P is the permutation matrix (same as in Dorig{l})
        nId = [nId1' Pmatrix3D(N,2)*nId1' Pmatrix3D(N,3)*nId1' Pmatrix3D(N,4)*nId1'];

end 