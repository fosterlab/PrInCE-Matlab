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



%% user.majorproteingroupsfile

fid = fopen(user.majorproteingroupsfile);

% Check csv format.
if ~strcmpi(user.majorproteingroupsfile(end-2:end),'csv')
  error('The following MQ file is not in .csv format:\n %s', user.majorproteingroupsfile)
end

% Check header
head = fgetl(fid);
head = strsplit(head,',');
if ~strcmpi(head{1},'majority protein ids')
  error('First line of the following file must be "majority protein IDs":\n %s', user.majorproteingroupsfile)
end

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



%% user.corumfile
if ~isempty(user.corumfile)
  
  fid = fopen(user.corumfile);
  
  % Check csv format.
  if ~strcmpi(user.corumfile(end-2:end),'csv')
    error('The following MQ file is not in .csv format:\n %s', user.corumfile)
  end
  
  % Check header
  head = fgetl(fid);
  head = strsplit(head,';');
  if length(head)==1
    warning('The following file should be ";" separated, but is likely "," separated:\n %s', user.majorproteingroupsfile)
    head = strsplit(head{1},',');
  end
  if isempty(strfind(lower(head{1}),'complex id'))
    error('First column of the following file must be "Complex id":\n %s', user.majorproteingroupsfile)
  end
  if isempty(strfind(lower(head{2}),'complex name'))
    error('Second column of the following file must be "Complex name":\n %s', user.majorproteingroupsfile)
  end
  if isempty(strfind(lower(head{5}),'uniprot'))
    error('Fifth column of the following file must contain "uniprot":\n %s', user.majorproteingroupsfile)
  end
  
  fclose(fid);
  
end

