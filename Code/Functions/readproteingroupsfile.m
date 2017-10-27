function Protein_IDs = readproteingroupsfile(user)
%READPROTEINGROUPSFILE Reads protein groups input file for PRINCE.
%   IDs=READFINALINTERACTIONLIST(user) reads the protein groups file
%   specified in the user structure. user is made by the master script
%   PRINCE.m. Protein groups file is a csv file where each row is a major
%   protein group, and each column in the row is a protein ID identified in
%   that group. IDs is a cell array with one protein ID per cell.
%
%   See also PRINCE.

% Very simple function.
% Reads user.majorproteingroupsfile and returns Protein_IDs.
% Protein_IDs is a cell array where each element is a protein name.
% Assumes user.majorproteingroupsfile is tab-delimited.
% Assumes user.majorproteingroupsfile has a header in the first line.

tmp = importdata(user.majorproteingroupsfile);
if isstruct(tmp)
  if isfield(tmp,'Sheet1')
    tmp2 = tmp.Sheet1;
  else
    tmp2 = tmp;
  end
elseif iscell(tmp)
  tmp2 = tmp;
end

Nprot = size(tmp2,1);

% Check if user.majorproteingroupsfile has a header
if ~isempty(strfind(tmp2{1},'rotein')) || ~isempty(strfind(tmp2{1},'roup'))
  istart = 2;
else
  istart = 1;
end

% Check if tmp2 already has one protein name per cell.
% To do this, check if there are any commas in any of the first cells.
Ncommas = zeros(Nprot,1);
for ii = istart:Nprot
  Ncommas(ii) = sum(double(tmp2{ii,1}==44));
end

% If there were no tabs or commas, tmp2 is correctly formatted.
% Clean it, assign it to Protein_IDs, and quit.
if sum(Ncommas)==0
  
  % Remove bad characters
  badchars = [',' '(' ')' '"'];
  for ii = 1:length(tmp2)
    I = ismember(tmp2{ii},badchars);
    tmp2{ii}(I) = [];
  end
  
  Protein_IDs = tmp2;

  
% If there were  tabs, tmp2 needs to be split up.
% Split each row of tmp2 according to tabs.
else
  
  Protein_IDs = cell(Nprot, 1000);
  Nmembers = 0;
  Nprot2 = 1;
  for ii = istart:Nprot
    
    % check if line has commas
    I = find(double(tmp2{ii,1}) == 44);
    
    if isempty(I)
      
      Protein_IDs{ii,1} = tmp2{ii,1};
    else
      
      consecutiveTabs = find(diff(I)==1)+1;
      
      Nprot2 = length(I) - length(consecutiveTabs);
      lastElem = I(end)-1;
      if ~isempty(consecutiveTabs)
        lastElem = I(consecutiveTabs(1)-1)-1;
      end
      
      if Nprot2 == 1
        Protein_IDs{ii,1} = tmp2{ii,1}(1:lastElem);
      elseif Nprot2 == 2
        Protein_IDs{ii,1} = tmp2{ii,1}(1:I(1)-1);
        Protein_IDs{ii,2} = tmp2{ii,1}(I(1)+1:lastElem);
      else
        
        Protein_IDs{ii,1} = tmp2{ii,1}(1:I(1)-1);
        for jj = 1:Nprot2 - 2
          Protein_IDs{ii,jj+1} = tmp2{ii,1}(I(jj)+1:I(jj+1)-1);
        end
        Protein_IDs{ii,Nprot2} = tmp2{ii,1}(I(Nprot2-1)+1:lastElem);
        
        % sanity check
        if jj+2 ~= Nprot2 && Nprot2~=2
          error('readproteingroupsfile: protein groups file improperly formatted.')
        end
      end
      
    end
    
    Nmembers = max([Nmembers Nprot2]);
  end
  
  Protein_IDs = Protein_IDs(:,1:Nmembers);
end
