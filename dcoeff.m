function [ res ] = dcoeff(coeff,s)


tmp = coeff;
for c=1:size(coeff,1)
    tmp(c,2)=coeff(c,s+2)*coeff(c,2);
    tmp(c,s+2)=coeff(c,s+2)-1;
end

res = [];
for c=1:size(tmp)
    if tmp(c,2)~=0
        res = cat(1,res,tmp(c,:));
    end
end