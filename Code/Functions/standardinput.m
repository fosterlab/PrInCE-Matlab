function standardinput(user)
% This function checks that input files specified by user have the correct format.
% - Check csv format.
% - Check header.
% - Check first 2 lines.
% - Check last 2 lines.

numeric_chars = [num2str(0:10) '-' '.'];

%% user.MQfiles
for ii = 1:length(user.MQfiles)
  fid = fopen(user.MQfiles{ii});
  
  % Check csv format.
  if ~strcmpi(user.MQfiles{ii}(end-2:end),'csv')
    error('The following MQ file is not in .csv format:\n %s', user.MQfiles{ii})
  end
  
  % Check header.
  head = fgetl(fid);
  head = strsplit(head,',');
  head = strtrim(head);
  head = head(~cellfun('isempty',head));
  if isempty(strfind(lower(head{1}),'protein'))
    fprintf('\nFirst column of the following MQ file must be protein IDs.')
    error('First row of the first column must contain the word "protein".\n %s', user.MQfiles{ii})
  end
  if isempty(strfind(lower(head{2}),'replicate'))
    fprintf('\nThe second column of the following MQ file must be the replicate number.')
    error('The first line of the second column must contain the word "replicate".\n %s', user.MQfiles{ii})
  end
  if length(head)-2 ~= user.Nfraction
    error('There is a mismatch between user.Nfraction (%d) and number of fractions in the MQ file (%d):\n %s',...
      user.Nfraction, length(head)-2, user.MQfiles{ii})
  end
  
  % Check first 2 lines.
  for jj = 1:2
    current_line = fgetl(fid);
    current_line = strsplit(current_line,',');
    current_line = strtrim(current_line);
    current_line = current_line(~cellfun('isempty',current_line));
    if length(current_line)-2 ~= user.Nfraction
      error('There is a mismatch between user.Nfraction (%d) and number of fractions in the MQ file (%d):\n %s',...
        user.Nfraction, length(head)-2, user.MQfiles{ii})
    end
    if ~ischar(current_line{1})
      error('First column of the following MQ file must be protein IDs.:\n %s',user.MQfiles{ii})
    end
    if ~isnumeric(current_line{2})
      if sum(ismember(current_line{2},numeric_chars))~=length(current_line{2})
        error('Second column of the following MQ file must be replicate number.:\n %s',user.MQfiles{ii})
      end
    end
  end
  
  % Check last 2 lines.
  current_line_m1 = [];
  while ~isequal(current_line, -1)
    current_line_m2 = current_line_m1;
    current_line_m1 = current_line;
    current_line = fgetl(fid);
  end
  clear current_line
  current_line{1} = strsplit(current_line_m2,',');
  current_line{2} = strsplit(current_line_m1,',');
  for jj = 1:2
    current_line{jj} = strtrim(current_line{jj});
    current_line{jj} = current_line{jj}(~cellfun('isempty',current_line{jj}));
    if length(current_line{jj})-2 ~= user.Nfraction
      error('There is a mismatch between user.Nfraction (%d) and number of fractions in the MQ file (%d):\n %s',...
        user.Nfraction, length(head)-2, user.MQfiles{ii})
    end
    if ~ischar(current_line{jj}{1})
      error('First column of the following MQ file must be protein IDs.:\n %s',user.MQfiles{ii})
    end
    if ~isnumeric(current_line{jj}{2})
      if sum(ismember(current_line{jj}{2},numeric_chars))~=length(current_line{jj}{2})
        error('Second column of the following MQ file must be replicate number.:\n %s',user.MQfiles{ii})
      end
    end
  end
  
  fclose(fid);
end

% ensure that user.Nreplicate is correct
tmp = readchromatogramfile2(user.MQfiles{1});
Nreplicate_empirical = length(unique(tmp.data(:,1)));
if user.Nreplicate ~= Nreplicate_empirical
  warning('Number of replicates is likely wrong in experimental_design.rtf.')
  warning(['Setting number of replicates to ' num2str(Nreplicate_empirical) '.'])
  user.Nreplicate = Nreplicate_empirical;
end
  


%% user.majorproteingroupsfile

if ~isempty(user.majorproteingroupsfile)
  
  fid = fopen(user.majorproteingroupsfile);
  
  % Check csv format.
  if ~strcmpi(user.majorproteingroupsfile(end-2:end),'csv')
    error('The following MQ file is not in .csv format:\n %s', user.majorproteingroupsfile)
  end
  
  % % Check header
  head = fgetl(fid);
  head = strsplit(head,',');
  % if ~strcmpi(head{1},'majority protein ids')
  %   error('First line of the following file must be "majority protein IDs":\n %s', user.majorproteingroupsfile)
  % end
  
  % Check all lines
  while ~feof(fid)
    current_line = fgetl(fid);
    current_line = strsplit(current_line,',');
    for jj = 1:length(current_line)
      if isempty(current_line{jj})
        continue;
      end
      if length(current_line{jj})>12
        error('The following file appears to contain badly formatted protein IDs:\n %s', user.majorproteingroupsfile)
      end
    end
  end
  
  fclose(fid);
  
end



%% user.corumfile
if ~isempty(user.corumfile)
  
  fid = fopen(user.corumfile);
  
  % Check header
  head = fgetl(fid);
  head1 = strsplit(head,';');
  head2 = strsplit(head,'\t');
  
  if length(head1)==1 & length(head2)>1
    if isempty(strfind((head2{1}),'ComplexID'))
      error('First column of the following file must be "ComplexID":\n %s', user.majorproteingroupsfile)
    end
    if isempty(strfind((head2{2}),'ComplexName'))
      error('Second column of the following file must be "ComplexName":\n %s', user.majorproteingroupsfile)
    end
    if isempty(strfind((head2{6}),'UniProt'))
      error('Sixth column of the following file must contain "UniProt":\n %s', user.majorproteingroupsfile)
    end
  elseif length(head1)>1 & length(head2)==1
    if isempty(strfind(lower(head1{1}),'complex id'))
      error('First column of the following file must be "Complex id":\n %s', user.majorproteingroupsfile)
    end
    if isempty(strfind(lower(head1{2}),'complex name'))
      error('Second column of the following file must be "Complex name":\n %s', user.majorproteingroupsfile)
    end
    if isempty(strfind(lower(head1{5}),'uniprot'))
      error('Fifth column of the following file must contain "uniprot":\n %s', user.majorproteingroupsfile)
    end
  end
  
  fclose(fid);
  
end

