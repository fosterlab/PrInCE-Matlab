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

% Remove self interactions
%M(eye(size(M,1))==1) = 0;

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
M2 = M>0;
cc = 0;
Icomps = find(sum(M2,2)>2);
clusters = cell(length(Icomps),1);
for uu = 1:length(Icomps)
    cc = cc+1;
    clusters{cc} = find(M2(Icomps(uu),:)>0);
end

% Merge equal clusters
tmp = nan(length(clusters));
for ii = 1:length(clusters)
    for jj = 1:length(clusters)
        if ii>=jj; continue; end
        tmp(ii,jj) = isequal(clusters{ii}, clusters{jj});
    end
end
for ii = 1:length(clusters)
    Iremove = find(tmp(ii,:) == 1);
    for jj = 1:length(Iremove)
    clusters{Iremove(jj)} = '-1';
    end
end
Iremove = cell2mat(cellfun(@(x) isequal(x,'-1'),clusters, 'UniformOutput', 0));
clusters(Iremove) = [];

