function user_new = standardinput(user)
% This function checks that input files specified by user have the correct format.
% - Check csv format.
% - Check header.
% - Check first 2 lines.
% - Check last 2 lines.

numeric_chars = [num2str(0:10) '-' '.'];

user_new = user;


%% user.MQfiles
cc = 0;
countedFraction = nan(10^5,1); % number of fractions in each row
replicate = nan(10^5,1); % replicate seen in each row
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
    warning('\nThe first column of the data file %s must be protein IDs.', user.MQfiles{ii})
  end
  if isempty(strfind(lower(head{2}),'replicate'))
    warning('\nThe second column of the data file %s must be the replicate number.', user.MQfiles{ii})
  end
  
  % Check first 2 lines.
  while not(feof(fid))
    cc = cc+1;
    current_line = fgetl(fid);
    current_line = strsplit(current_line,',');
    current_line = strtrim(current_line);
    current_line = current_line(~cellfun('isempty',current_line));
    countedFraction(cc) = length(current_line) - 2;
    replicate(cc) = str2double(current_line{2});
    if ~ischar(current_line{1})
      error('First column of the following MQ file must be protein IDs.:\n %s',user.MQfiles{ii})
    end
    if ~isnumeric(current_line{2})
      if sum(ismember(current_line{2},numeric_chars))~=length(current_line{2})
        error('Second column of the following MQ file must be replicate number.:\n %s',user.MQfiles{ii})
      end
    end
  end
  
    
  fclose(fid);
end

% Count how many fractions were detected
countedFraction = countedFraction(1:cc);
fractionCounts = unique(countedFraction);
goodFractions = 0;
if length(fractionCounts)==1
    % If all rows had the same number of fractions
    user_new.Nfraction = fractionCounts;
    goodFractions = 1;
else
    % If there were different fraction numbers
    user_new.Nfraction = max(fractionCounts);
    ss = sprintf('The number of FRACTIONS is not consistent in condition files. \nUsing the max number of fractions (%s). \nIf %s is incorrect, check that files are correctly formatted.',...
        num2str(user_new.Nfraction),num2str(user_new.Nfraction));
    warning(ss);
    hh = warndlg(sprintf(ss));
    %uiwait(hh)
end

% Count how many replicates were detected
replicate = replicate(1:cc);
replicate = unique(replicate);
replicate(isnan(replicate)) = [];
replicate(isempty(replicate)) = [];
user_new.Nreplicate = length(replicate);
goodReplicates = 0;
if user_new.Nreplicate>15
    % If too many replicates were detected
    ss = sprintf('%s REPLICATES detected in %s. Thats a lot! \nCheck file to ensure replicate column is correctly formatted.',...
        num2str(user_new.Nreplicate), [user_new.silacratios{ii} '.csv']);
    warning(ss);
    hh = warndlg(sprintf(ss));
    uiwait(hh)
else
    % If a reasonable number were detected
    goodReplicates = 1;
end

% Tell the user how many replicates/fractions were detected
if goodReplicates==1 && goodFractions ==1
    ss = sprintf('Detected %s FRACTIONS and %s REPLICATES in condition files. If that is not correct, check that condition files are correctly formatted.',...
        num2str(user_new.Nfraction), num2str(user_new.Nreplicate));
    sprintf(ss)
    msgbox(ss, 'Replicates and fractions');
elseif goodReplicates==1 && goodFractions == 0
    ss = sprintf('Detected %s REPLICATES in condition files. If that is not correct, check that condition files are correctly formatted.',...
        num2str(user_new.Nreplicate));
    sprintf(ss)
    msgbox(ss, 'Replicates');
elseif goodReplicates==0 && goodFractions == 1
    ss = sprintf('Detected %s FRACTIONS in condition files. If that is not correct, check that condition files are correctly formatted.',...
        num2str(user_new.Nfraction));
    sprintf(ss)
    msgbox(ss, 'Fractions');
end


%% user.majorproteingroupsfile

if ~isempty(user.majorproteingroupsfile)
  fid = fopen(user.majorproteingroupsfile);
  
  % Check csv format.
  if ~strcmpi(user.majorproteingroupsfile(end-2:end),'csv')
    error('The following MQ file is not in .csv format:\n %s', user.majorproteingroupsfile)
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
        warning('The following file appears to contain badly formatted protein IDs:\n %s', user.majorproteingroupsfile)
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
      warning('First column of the following file should be "ComplexID":\n %s', user.majorproteingroupsfile)
    end
    if isempty(strfind((head2{2}),'ComplexName'))
      warning('Second column of the following file should be "ComplexName":\n %s', user.majorproteingroupsfile)
    end
    if isempty(strfind((head2{6}),'UniProt'))
      warning('Sixth column of the following file should contain "UniProt":\n %s', user.majorproteingroupsfile)
    end
  elseif length(head1)>1 & length(head2)==1
    if isempty(strfind(lower(head1{1}),'complex id'))
      warning('First column of the following file should be "Complex id":\n %s', user.majorproteingroupsfile)
    end
    if isempty(strfind(lower(head1{2}),'complex name'))
      warning('Second column of the following file should be "Complex name":\n %s', user.majorproteingroupsfile)
    end
    if isempty(strfind(lower(head1{5}),'uniprot'))
      warning('Fifth column of the following file should contain "uniprot":\n %s', user.majorproteingroupsfile)
    end
  end
  
  fclose(fid);
end

