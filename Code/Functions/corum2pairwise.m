function corum2pairwise(user)

%CORUM2PAIRWISE Creates a csv list of pairwise interactions
%    using the CORUM file allComplexes.csv.
%
%    Original script by Anders Kristensen.
%    Modified by Nichollas Scott, 2014.
%    Modified by Greg Stacey, 2016.



%% Read raw CORUM file
fid = fopen(user.corumfile,'r');

% check header is correctly formatted
header = fgetl(fid);
if isempty(strfind(header,'Complex id')) || isempty(strfind(header,'subunits (UniProt IDs)'))
  error('corum2pairwise: Improperly formatted corum file (missing header?).')
end

% read body of file
cc = 0;
organism = cell(10000,1);
complexes = cell(10000,1);
while 1
  t = fgetl(fid);
  if(~ischar(t)),break,end
  cc = cc+1;
  Idelim = strfind(t, ';');
  
%   Iquotes = strfind(t,'"');
%   % remove ',' between each pair of quotes
%   for ii = 1:length(Iquotes)/2
%     i1 = Iquotes((ii-1)*2 +1);
%     i2 = Iquotes(ii*2);
%     Idelim(Idelim > i1 & Idelim < i2) = [];
%   end
  
  organism{cc} = t(Idelim(3)+1 : Idelim(4)-1);
  complexes{cc} = t(Idelim(4)+1 : Idelim(5)-1);
  complexes{cc} = strrep(complexes{cc}, '"', '');
end
organism = organism(1:cc);
complexes = complexes(1:cc);

fclose(fid);



%% Split complexes into pairwise list

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
  Idelim = [0 strfind(cmplx, ',') length(cmplx)+1];
  Nprot = length(Idelim) - 1;
  if Nprot<2
    continue
  end
  
  for jj = 1:Nprot
    prot1 = cmplx(Idelim(jj)+1 : Idelim(jj+1)-1);
    for kk = 1:Nprot
      if jj>=kk % only scan the upper-triangular matrix
        continue
      end
      
      cxi = cxi+1;
      prot2 = cmplx(Idelim(kk)+1 : Idelim(kk+1)-1);
      
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



%% Write pairwise interactions to file

fn = [user.maindir '/Data/Corum_pairwise.csv'];
fid = fopen(fn,'w');
for ii = 1:size(pairwiselist,1)
  fprintf(fid,'%s,%s\n',pairwiselist{ii,1},pairwiselist{ii,2});
end
fclose(fid);


