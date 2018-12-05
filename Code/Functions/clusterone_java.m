function Members = clusterone_java(M, pp, density_threshold, java_path)

%CLUSTERONE_JAVA Creates a list of protein complexes using the
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

this_path = pwd;

% parameters
if nargin<2
  density_threshold = 2;
  pp = 0;
end

% Sanity checks
if size(M,1)~=size(M,2)
  error('myclusterone: Interaction matrix must be square')
end
if sum(M(:)<0) >= 1
  error('myclusterone: Interaction matrix must contain only non-negative entries')
end


% 1. write M to file
[ia,ib] = find(M>0);
fn_tmpM = [this_path '/.Mtmp.txt'];
fid = fopen(fn_tmpM, 'w');
for ii = 1:length(ia)
  fprintf(fid,'%d %d %6.4f\n',ia,ib,M(ia(ii),ib(ii)));
end
fclose(fid);


% 2. Call clusterone java
fn_tmpout = [this_path '/.outtmp.txt'];
system_call = ['java -jar ' java_path ' ' fn_tmpM ' -s 2 -d ' num2str(density_threshold) ...
  ' --penalty ' num2str(pp) ' > ' fn_tmpout];
system(system_call);


% 3. Read clusterone java output
fid = fopen(fn_tmpout);
Members = cell(10^5,1);
cc = 0;
while ~feof(fid)
  cc = cc+1;
  t1 = strsplit(fgetl(fid),'\t');
  Members{cc} = nan(length(t1));
  for ii = 1:length(t1)
    Members{cc}(ii) = str2double(t1{ii});
  end
end
Members = Members(1:cc);

% 4. Clean up
delete(fn_tmpM) 
delete(fn_tmpout)


