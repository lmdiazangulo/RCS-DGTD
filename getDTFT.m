function [S] = getDTFT(fq, t, s)

N = numel(s);

S = conj(exp(- 1i *2*pi*fq'*t)*s' / N);