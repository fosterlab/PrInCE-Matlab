function [Members, Density] = myclusterone(M, pp, density_threshold, iterMax)

%MYCLUSTERONE Creates a list of protein complexes using the
%    ClusterONE algorithm (Nepusz 2012).
%    [Members,Connections] = MYCLUSTERONE(M,P) groups the 
%    interactions in pairwise interaction matrix M. Members is 
%    an Nx1 cell array, where each cell contains the indices 
%    of all proteins in a single complex. Connections is an 
%    a square matrix, of which each non-zero entry denotes a
%    interaction. The proteins associated with each row/column 
%    of Connections matrix are given by the corresponding 
%    Members entry.
%
%    Three parts:
%    1. greedy growth
%    2. merge overlapping groups
%    3. discard groups with <3 members or low density
%
% Greedy growth algorithm:
% while ~all_proteins_are_in_a_complex
%   V = seed = protein w/ most vertices not already in a complex
%   t=0
%   while grow
%     calculate cohesivness of Vt, f(Vt)
%     V_{t+1} = Vt
%     make ve, which is "every external vertex incident on at least one boundary edge"
%     if f(Vt union ve) > f(Vt), V_{t+1} = Vt union ve
%     make vi, which is "ever internal vertex incident on at least one boundary edge"
%     if f(Vt / vi) > f(Vt), V_{t+1} = Vt / vi
%     grow = Vt == V_{t+1}
%
% ClusterONE ref:
% Detecting overlapping protein complexes in protein-protein 
% interaction networks Nature Methods 9, 471?472 (2012) 
% doi:10.1038/nmeth.1938 


% parameters
if nargin<2
  density_threshold = 2;
  pp = 0;
end

if nargin<4
    iterMax = 50;
end

% Sanity checks
if size(M,1)~=size(M,2)
  error('myclusterone: Interaction matrix must be square')
end
if sum(M(:)<0) >= 1
  error('myclusterone: Interaction matrix must contain only non-negative entries')
end

  
% Count how many interactions each node is in
Nprot = size(M,1);
Nconnections = sum(M>0,2);

% Binary vector
% Used to keep track of which proteins have already been assigned to a complex
inacomplex = zeros(Nprot,1);

% 1. Grow complexes
Members = cell(100000,1);
cmplxcount = 0;
while sum(inacomplex)<Nprot && cmplxcount<size(M,1)

  % Choose protein with most number of connections as starting seed, V0
  Nconnections(inacomplex==1) = -1;
  [~,V] = max(Nconnections);
  inacomplex(V) = 1;
  
  grow = 1;
  kk = 0;
  while grow && kk<iterMax
    kk = kk+1;
    V0 = V;
    
    % V, current complex
    % find every boundary edge
    [vi0, ve0] = find(M(V,:)>0);
    vi0 = V(vi0);
    Inotboundary = ismember(ve0,V);
    vi0 = vi0(~Inotboundary); % the internal vertex of each boundary edge
    ve0 = ve0(~Inotboundary); % the external vertex of each boundary edge
    coh0 = sum(sum(M(V,V))) / (sum(sum(M(V,V))) + sum(sum(M(vi0,ve0))) + pp*length(V));
    if size(ve0,1)>size(ve0,2)
      ve0 = ve0';
    end
    if size(vi0,1)>size(vi0,2)
      vi0 = vi0';
    end
    
    % V1, V0 + external vertices
    V1 = unique([V ve0]);
    [vi, ve] = find(M(V1,:)>0);
    vi = V1(vi);
    Inotboundary = ismember(ve,V1);
    vi = vi(~Inotboundary); % the internal vertex of each boundary edge
    ve = ve(~Inotboundary); % the external vertex of each boundary edge
    coh1 = sum(sum(M(V1,V1))) / (sum(sum(M(V1,V1))) + sum(sum(M(vi,ve))) + pp*length(V1));
    
    % V1, V0 - internal vertices
    V2 = V(~ismember(V,vi0));
    [vi, ve] = find(M(V2,:)>0);
    vi = V2(vi);
    Inotboundary = ismember(ve,V2);
    vi = vi(~Inotboundary); % the internal vertex of each boundary edge
    ve = ve(~Inotboundary); % the external vertex of each boundary edge
    coh2 = sum(sum(M(V2,V2))) / (sum(sum(M(V2,V2))) + sum(sum(M(vi,ve))) + pp*length(V2));
    
    if coh1>coh0
      V = V1;
    end
    if coh2>coh0
      V = V2;
    end
    
    grow = ~isequal(V,V0);
  end
  
  cmplxcount = cmplxcount+1;
  Members{cmplxcount} = V;
  
  inacomplex(V) = 1;
end
Members = Members(1:cmplxcount);
% Remove complexes of length 1 and 2 
% These can't be merged using the threshold score of 0.8 and they slow things down.
for ii = 1:size(Members)
  if length( Members{ii})<2
    Members{ii} = [];
  end
end
Members = Members(~cellfun('isempty',Members));
cmplxcount = length(Members);

% 2. Merge complexes
% Calculate every complex-pairwise overlap score
overlapscore = zeros(cmplxcount,cmplxcount);
for ii = 1:cmplxcount
  for jj = 1:cmplxcount
    if ii==jj;continue;end
    overlapscore(ii,jj) = length(intersect(Members{ii},Members{jj}))^2 / length(Members{ii}) / length(Members{jj});
  end
end

overlapscore = triu(overlapscore);

[Ix, Iy] = find(overlapscore>0.8);
for ii = 1:length(Ix)
  Members{Ix(ii)} = unique([Members{Ix(ii)} Members{Iy(ii)}]);
  Members{Iy(ii)} = [];
  
  Ix(Ix==Iy(ii)) = Ix(ii);
  Iy(Ix==Iy(ii)) = Ix(ii);
end
Members = Members(~cellfun('isempty',Members));


% 3. Remove small or low-density complexes
Density = nan(size(Members));
for ii = 1:size(Members)
  I = Members{ii};
  m = M(I,I);
  n = length(I);
  Density(ii) = sum(m(:)) / n / n;
end
I = Density(ii)<density_threshold;
Density(I) = [];
Members(I) = [];
Members = Members(~cellfun('isempty',Members));
