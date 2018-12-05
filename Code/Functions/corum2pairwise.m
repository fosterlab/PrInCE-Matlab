function corum2pairwise(user)

%CORUM2PAIRWISE Split reference database into pairwise interactions.
%   CORUM2PAIRWISE(user) creates a csv file with two columns, using the 
%   complexes contained in file user.corumfile. Each complex of size N is 
%   split into N(N-1)/2 pairwise interactions. Created csv file is 
%   'Corum_pairwise.csv' and is located in the Output/tmp/ folder.
%
%   user is a structure and must be created by master script PRINCE, or 
%   otherwise have a similar format.
%
%   user.corumfile must be in the format of CORUM's allComplexes.txt or 
%   allComplexes.csv.
%
%   See also PRINCE




% Read raw CORUM file
fid = fopen(user.corumfile,'r');

% check header is correctly formatted
header = fgetl(fid);
corumType = 1;
if isempty(strfind(header,'Complex id')) || isempty(strfind(header,'subunits (UniProt IDs)'))
  if isempty(strfind(header,'ComplexName'))
    error('corum2pairwise: Improperly formatted corum file (missing header?).')
  else
    corumType = 2;
  end
end

% read body of file
if corumType == 1
  cc = 0;
  organism = cell(10000,1);
  complexes = cell(10000,1);
  while 1
    t = fgetl(fid);
    if(~ischar(t)),break,end
    cc = cc+1;
    Idelim = strfind(t, ';');
    
    organism{cc} = t(Idelim(3)+1 : Idelim(4)-1);
    complexes{cc} = t(Idelim(4)+1 : Idelim(5)-1);
    complexes{cc} = strrep(complexes{cc}, '"', '');
  end
  organism = organism(1:cc);
  complexes = complexes(1:cc);
  
elseif corumType == 2
  cc = 0;
  organism = cell(10000,1);
  complexes = cell(10000,1);
  while ~feof(fid)
    t1 = strsplit(fgetl(fid),'\t');
    cc = cc+1;
    
    organism{cc} = t1{3};
    complexes{cc} = t1{6};
  end
  organism = organism(1:cc);
  complexes = complexes(1:cc);
  
end

fclose(fid);



% Split complexes into pairwise list

pairwiselist = cell(10^6,2);
cxi = 0;
for ii = 1:cc
  % Check that the current complex is from the desired organism
  if isfield(user,'organism')
    if ~isempty(user.organism)
      if ~strcmpi(user.organism, organism{ii})
        
        % Organism can be 'Mammalia', which includes 'mouse', 'rat', and 'human'
        if strcmpi(organism{ii},'mammalia') && ...
            (strcmpi(user.organism,'mouse') || strcmpi(user.organism,'rat') || ...
            strcmpi(user.organism,'human') || strcmpi(user.organism,'pig'))
        else
          continue
        end
      end
    end
  end
  
  cmplx = complexes{ii};
  cmplx = strrep(cmplx,'"',' ');
  cmplx = strrep(cmplx,',',' ');
  cmplx = strrep(cmplx,';',' ');
  cmplx = strsplit(cmplx,' ');
  cmplx = cmplx(not(cellfun('isempty',cmplx)));
  cmplx = unique(cmplx);
  Nprot = length(cmplx);

  if Nprot<2
    continue
  end
  
  for jj = 1:Nprot
    prot1 = cmplx{jj};
    for kk = 1:Nprot
      if jj>=kk % only scan the upper-triangular matrix
        continue
      end
      
      cxi = cxi+1;
      prot2 = cmplx{kk};
      
      % remove parentheses from protein IDs
      prot1 = strrep(prot1,'(','');
      prot1 = strrep(prot1,')','');
      prot2 = strrep(prot2,'(','');
      prot2 = strrep(prot2,')','');
      
      pairwiselist{cxi,1} = prot1;
      pairwiselist{cxi,2} = prot2;
    end
  end
end
pairwiselist = pairwiselist(1:cxi,:);



% Write pairwise interactions to file

if ~exist([user.maindir '/Output'], 'dir'); mkdir([user.maindir '/Output']); end
if ~exist([user.maindir '/Output/tmp'], 'dir'); mkdir([user.maindir '/Output/tmp']); end
fn = [user.maindir '/Output/tmp/Corum_pairwise.csv'];
fid = fopen(fn,'w');
for ii = 1:size(pairwiselist,1)
  fprintf(fid,'%s,%s\n',pairwiselist{ii,1},pairwiselist{ii,2});
end
fclose(fid);


