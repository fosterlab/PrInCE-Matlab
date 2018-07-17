function M = mymcl(M,p,e,minval,iterMax)

if nargin<2
    p = 2;
end
if nargin<3
    e = 2;
end
if nargin<4
    minval = 0.001;
end
if nargin<5
    iterMax = 20;
end

ii = 0;
while ii<iterMax
    % normalize columns
    M = bsxfun(@rdivide, M, sum(M));
    
    % Expand by e
    M = mpower(M, e);
    
    % Inflate by p
    M = power(M, p);
    
    % normalize columns
    M = bsxfun(@rdivide, M, sum(M));
    
    % prune
    M(M<minval) = 0;
    
    ii = ii+1; 
end
