function user = standardinput(user)
% This function checks that input files specified by user have the correct format.
% - Check csv format.
% - Check header.
% - Check first 2 lines.
% - Check last 2 lines.

numeric_chars = [num2str(0:10) '-' '.'];

%user = user;


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
  user.Nfraction = fractionCounts;
  goodFractions = 1;
else
  % If there were different fraction numbers
  user.Nfraction = max(fractionCounts);
  ss = sprintf('The number of FRACTIONS is not consistent in condition files. \nUsing the max number of fractions (%s). \nIf %s is incorrect, check that files are correctly formatted.',...
    num2str(user.Nfraction),num2str(user.Nfraction));
  warning(ss);
  hh = warndlg(sprintf(ss));
  %uiwait(hh)
end

% Count how many replicates were detected
replicate = replicate(1:cc);
replicate = unique(replicate);
replicate(isnan(replicate)) = [];
replicate(isempty(replicate)) = [];
user.Nreplicate = length(replicate);
goodReplicates = 0;
if user.Nreplicate>15
  % If too many replicates were detected
  ss = sprintf('%s REPLICATES detected in %s. Thats a lot! \nCheck file to ensure replicate column is correctly formatted.',...
    num2str(user.Nreplicate), [user.silacratios{ii} '.csv']);
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
    num2str(user.Nfraction), num2str(user.Nreplicate));
  sprintf(ss)
  msgbox(ss, 'Replicates and fractions');
elseif goodReplicates==1 && goodFractions == 0
  ss = sprintf('Detected %s REPLICATES in condition files. If that is not correct, check that condition files are correctly formatted.',...
    num2str(user.Nreplicate));
  sprintf(ss)
  msgbox(ss, 'Replicates');
elseif goodReplicates==0 && goodFractions == 1
  ss = sprintf('Detected %s FRACTIONS in condition files. If that is not correct, check that condition files are correctly formatted.',...
    num2str(user.Nfraction));
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


%% user.MQfiles - are they identically ordered?
% This needs to done after setting user.Nfraction

% Check that condition files have identical proteins in the same order
% For now, read in entire protein ID list
clear chromdata
for ii = 1:length(user.MQfiles)
  chromdata(ii) = readchromatogramfile2(user.MQfiles{ii}, user.Nfraction);
  chromdata(ii).textdata = chromdata(ii).textdata(2:end,:);
end

% Check that condition files have identical proteins in the same order
protPairIdentical = nan(length(user.MQfiles));
for ii = 1:length(user.MQfiles)
  for jj = 1:length(user.MQfiles)
    if ii>=jj; continue; end
    protPairIdentical(ii,jj) = isequal(chromdata(ii).textdata, chromdata(jj).textdata);
  end
end
if nansum(protPairIdentical(:)) < (length(user.MQfiles) * (length(user.MQfiles)-1)) / 2
  disp('Writing temporary condition files.')
  
  % Condition files are not identically ordered by protein
  % Solution: Make temporary condition files with identical protein order
  
  % find all unique proteins
  for ii = 1:length(user.MQfiles)
    if ii==1
      allProts = unique(chromdata(ii).textdata);
    else
      allProts = unique([allProts; chromdata(ii).textdata]);
    end
  end
  
  % find all unique replicates
  for ii = 1:length(user.MQfiles)
    if ii==1
      allReps = unique(chromdata(ii).replicate);
    else
      allReps = unique([allReps; chromdata(ii).replicate]);
    end
  end
  
  % Write condition files
  for ii = 1:length(user.MQfiles)
    % make data
    nrows = length(allReps) * length(allProts);
    data2write.data = nan(nrows, user.Nfraction+1);
    data2write.text = cell(nrows,1);
    cc = 0;
    for jj = 1:length(allReps)
      rep = allReps(jj);
      I1 = chromdata(ii).replicate == rep;
      for kk = 1:length(allProts)
        cc = cc+1;
        prot = allProts{kk};
        data2write.text{cc} = prot;
        data2write.data(cc,1) = rep;
        
        % get chromatogram from original data
        I = ismember(chromdata(ii).textdata, prot) & I1;
        if sum(I)==0; continue; end
        if sum(I)~=1; error(['Duplicate protein IDs in the same replicate: ' user.MQfiles{ii}]); end
        
        data2write.data(cc,:) = [chromdata(ii).replicate(I) chromdata(ii).data(I,:)];
      end
    end
    
    % get original header
    fid = fopen(user.MQfiles{ii});
    head = fgetl(fid);
    fclose(fid);
    
    % write output
    fn = [user.maindir '/Output/tmp/tmp_' user.silacratios{ii} '.csv'];
    fid = fopen(fn,'w');
    fprintf(fid,'%s\n',head);
    for jj = 1:size(data2write.data,1)
      fprintf(fid,'%s,%d,', data2write.text{jj}, data2write.data(jj,1));
      for kk = 1:size(data2write.data,2)-1
        fprintf(fid,'%6.4f,', data2write.data(jj,kk+1));
      end
      fprintf(fid,'\n');
    end
    
    % rename user.MQfiles{ii}
    user.MQfiles{ii} = fn;
  end
end


% Are there duplicate protein IDs in the same replicate
for ii = 1:length(user.MQfiles)
  reps = unique(chromdata(ii).replicate);
  for jj = 1:length(reps)
    I = chromdata(ii).replicate==reps(jj);
    if length(unique(chromdata(ii).textdata)) < sum(I)
      error(['Duplicate protein IDs in the same replicate: ' user.MQfiles{ii}])
    end
  end
end


