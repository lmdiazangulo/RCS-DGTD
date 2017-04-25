function [res]=comb(m,n)
% function [res]=comb(m,n)
% discrete combinatorial number. m over n

res = factorial(m)/(factorial(n)*factorial(m-n));

end