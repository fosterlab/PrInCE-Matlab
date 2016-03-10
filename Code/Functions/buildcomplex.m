function [Members, Connections] = buildcomplex(M)

%BUILDCOMPLEX Creates a list of protein complexes.
%    [Members,Connections] = pcpsilacmcl(M,P) groups the 
%    interactions in pairwise interaction matrix M. Members is 
%    an Nx1 cell array, where each cell contains the indices 
%    of all proteins in a single complex. Connections is an 
%    a square matrix, of which each non-zero entry denotes an
%    interaction. The proteins associated with each row/column 
%    of Connections matrix are given by the corresponding 
%    Members entry.
%    
%    Algorithm: 
%    *for each protein A
%      *find all proteins (B) directly interacting
%      *for each protein B
%        *note that A and B are in a complex
%        *find all proteins (C) directly interacting with B
%        *for each protein C, repeat
%        *stop when no new proteins are found
%    Find all direct interactions with protein A. Find all 
%    direct interactions with those proteins. Continue
%    iterating like this until no new proteins found. Call
%    that group a complex. Procede to the next protein not
%    already in a complex.
%
% Created by Greg Stacey, Feb 2016.


Nprot = size(M,1);

% Reduce M to just upper-triangular
M = triu(M);

% Binary vector
% Used to keep track of which proteins have already been assigned to a complex
bv = zeros(Nprot,1);

% Start finding complexes
Members = cell(1000,1);
Connections = cell(1000,1);
countcmplx = 0;
for jj = 1:Nprot
  % check if this protein is already in a complex
  if bv(jj)==1
    continue;
  end
  
  %I = find(M(jj,:)>=minrep | M(:,jj)'>=minrep);
  I = find(M(jj,:)>0 | M(:,jj)'>0);
  if ~isempty(I)
    % initialize this complex
    countcmplx = countcmplx+1;
    bv(jj) = 1;
    membs = [jj I];
    
    % BOOM! Another option is to IGNORE open_branches. Just keep exploring ALL branches, adding
    % interactions to the list, and pruning with unique. When the list stops growing, you've
    % stopped finding new interactions, and the complex is complete. Ignore the binary vector!
    
    % start exploring connections
    open_branches = membs(bv(membs)==0);
    %delta_membs = 1;
    while ~isempty(open_branches)
      %while delta_membs>0
      %nmembs = length(membs);
      for kk = 1:length(open_branches)
        I2 = find(M(open_branches(kk),:)>0 | M(:,open_branches(kk))'>0);
        membs = [membs I2];
      end
      membs = unique(membs);
      open_branches = membs(bv(membs)==0);
      bv(open_branches) = 1;
      %delta_membs = length(membs) - nmembs;
      %[length(open_branches) delta_membs]
    end
    
    Members{countcmplx} = membs;
    Connections{countcmplx} = M(membs,membs);
  end
end

Members = Members(1:countcmplx);
Connections = Connections(1:countcmplx);
