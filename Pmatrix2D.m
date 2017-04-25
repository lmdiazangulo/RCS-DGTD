function [ P ] = Pmatrix2D(N)
% P = Rotation Matrix

Np = (N+1)*(N+2)/2;
Nfp = N+1;

originalNum = 1;
original = zeros(Nfp,Nfp);
for i = 1:Nfp
    for j = 1:i
        original(i,j) = originalNum;
        originalNum = originalNum + 1;
    end
end
 
rotatedNum = 1;
rotated = zeros(Nfp,Nfp);
for i = 1:Nfp
    j = Nfp-1;
    while j >= Nfp - i % This can be done with a for with A:(-1):B
        rotated(j+1, (j-Nfp+i+1)) = rotatedNum;
        rotatedNum = rotatedNum + 1;
        j = j-1;
    end
end

P = zeros(Np,Np);
for i = 1:Nfp
    for j = 1:i
        P( rotated(i,j),  original(i,j)) = 1;
    end
end