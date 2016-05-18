function corumextractcomplex(user)

%CORUMEXTRACTCOMPLEX Creates a csv list of complexes
%    using the CORUM file allComplexes.csv.



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
  if sum(Idelim)==0
    Idelim = strfind(t, ',');
  end
  organism{cc} = t(Idelim(3)+1 : Idelim(4)-1);
  complexes{cc} = t(Idelim(4)+1 : Idelim(5)-1);
end
organism = organism(1:cc);
complexes = complexes(1:cc);

fclose(fid);



%% Remove complexes from irrelevant organisms

for ii = 1:cc
  
  % Check that the current complex is from the desired organism
  if isfield(user,'organism')
    if ~isempty(user.organism)
      if ~strcmpi(user.organism, organism{ii})
        
        % Organism can be 'Mammalia', which includes 'mouse', 'rat', 'human', and 'pig
        if strcmpi(organism{ii},'mammalia') && ...
            (strcmpi(user.organism,'mouse') || strcmpi(user.organism,'rat') || ...
            strcmpi(user.organism,'human') || strcmpi(user.organism,'pig'))
        else
          complexes{ii} = [];
        end
      end
    end
  end
  
  complexsize = length(strfind(complexes{ii}, ','));
  if complexsize<1
    complexes{ii} = [];
  end
end
I = ~cellfun('isempty',complexes);
complexes = complexes(I);



%% Write complexes to file

fn = [user.maindir '/Output/tmp/Corum_complex.csv'];
fid = fopen(fn,'w');
for ii = 1:size(complexes,1)
  cmplx = complexes{ii};
  
  Idelim = [0 strfind(cmplx, ',') length(cmplx)+1];
  Nprot = length(Idelim) - 1;
  if Nprot<2
    error('corumextractcomplex: Complex of size 1 detected.')
  end
  
  for jj = 1:Nprot-1
    prot1 = cmplx(Idelim(jj)+1 : Idelim(jj+1)-1);
    
    % remove parentheses from protein IDs
    prot1 = strrep(prot1,'(','');
    prot1 = strrep(prot1,')','');
    
    fprintf(fid,'%s,',prot1);
  end
  jj = Nprot;
  prot1 = cmplx(Idelim(jj)+1 : Idelim(jj+1)-1);
  
  % remove parentheses from protein IDs
  prot1 = strrep(prot1,'(','');
  prot1 = strrep(prot1,')','');
  
  fprintf(fid,'%s',prot1);
  fprintf(fid,'\n');
end
fclose(fid);


