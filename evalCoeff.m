function [ res ] = evalCoeff(coeff,p)

numeval = size(p,1);
dim = size(p,2);

res = zeros(numeval,1);
for j=1:numeval
    switch (dim)
        case 3
            for i=1:size(coeff,1);
                res(j)=res(j)+coeff(i,2)* p(j,1)^coeff(i,3)* p(j,2)^coeff(i,4)* p(j,3)^coeff(i,5);
            end
        otherwise
            disp('Problem with evalCoeff');
    end
end
        