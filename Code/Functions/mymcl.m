function [M, clusters] = mymcl(M,p,e,minval,iterMax)

% Takes
%   connection matrix M
%   parameters p, e, minval, iterMax.
% Returns
%   MCL-matrix
%   clusters in a cell-array

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


% Reduce M to remove rows/columns with zero sum
Inotzero = find(not(sum(M)==0));
Mnozeros = M(Inotzero,Inotzero);


ii = 0;
while ii<iterMax
    % normalize columns
    Mnozeros = bsxfun(@rdivide, Mnozeros, sum(Mnozeros));
    
    % Expand by e
    Mnozeros = mpower(Mnozeros, e);
    
    % Inflate by p
    Mnozeros = power(Mnozeros, p);
    
    % normalize columns
    Mnozeros = bsxfun(@rdivide, Mnozeros, sum(Mnozeros));
    
    % prune
    Mnozeros(Mnozeros<minval) = 0;
    
    ii = ii+1;
end


% Add zero-sum rows/columns back in
M = zeros(size(M));
M(Inotzero, Inotzero) = Mnozeros;


% Turn MCL-matrix M into a cell array of clusters
M = M>0;
cc = 0;
Icomps = find(sum(M,2)>2);
clusters = cell(length(Icomps),1);
for uu = 1:length(Icomps)
    cc = cc+1;
    clusters{cc} = find(M(Icomps(uu),:)>0);
end
