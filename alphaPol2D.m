function [ alpha ] = alphaPol2D(r1, r2, r3)

alpha = zeros(numel(r1),numel(r2),numel(r3));
for i=1:numel(r2)
    for j=1:numel(r3)
        alpha(:,i,j) = conv( conv(r1,r2(i)), r3(j));
    end
end

end
