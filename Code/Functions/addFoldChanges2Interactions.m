% Goal: Write these columns in Final_interactions.csv:
%   1. p-value (fold changes over replicate)
%   2. avg fold change
%   3. fold change in each replicate
%
% Workflow:
%   1. Find fold changes w/in replicate
%   2. Align b/w replicates
%   3. Combine fold changes across replicates
%   4. Collate fold changes for a single protein
%   5. Write fold changes to final_interactions.csv




%% 0. Initialize

maindir = user.maindir;

% get raw chromatograms
rawdatafiles = dir([maindir 'Input/condition*']);
clear rawdata
for ii = 1:length(rawdatafiles)
    ii
    tmp = readchromatogramfile2([maindir 'Input/' rawdatafiles(ii).name]);
    % remove 'sheet1' fields
    if isfield(tmp,'Sheet1')
        tmp = tmp.Sheet1;
    end
    fn3 = fieldnames(tmp);
    for jj = 1:length(fn3)
        tmp1 = tmp.(fn3{jj});
        if isfield(tmp1,'Sheet1')
            tmp.(fn3{jj}) = tmp.(fn3{jj}).Sheet1;
        end
    end
    
    rawdata(ii).data = tmp.data;
    rawdata(ii).Proteins = tmp.textdata;
    
    % if txt_val includes protein groups, reduce it to the first protein in each group
    for jj = 1:size(rawdata(ii).Proteins,1)
        tmp = strsplit(rawdata(ii).Proteins{jj},';');
        rawdata(ii).Proteins{jj} = tmp{1};
    end
    
    % if rawdata & txt_val are the same length, assume they both have headers, remove rawdata header
    if size(rawdata(ii).data,1)==size(rawdata(ii).Proteins,1)
        rawdata(ii).data = rawdata(ii).data(2:end,:);
    end
end


%%   1. Find fold changes w/in replicate
% This is done by FoldChanges.m

% Read fold changes
fn1 = [maindir 'Output/Data/FoldChanges/Fold_change_condition1_vs_condition2.csv'];
fid = fopen(fn1);
head = fgetl(fid); % skip header
clear FC
FC.protein = cell(10^5,1);
FC.replicate = nan(10^5,1);
FC.fraction = nan(10^5,1);
FC.foldchange = nan(10^5,1);
cc = 0;
while not(feof(fid))
    t1 = strsplit(fgetl(fid),',','collapsedelimiters',0);
    cc = cc+1;
    FC.protein{cc} = t1{1};
    FC.replicate(cc) = str2double(t1{2});
    FC.fraction(cc) = str2double(t1{3});
    %FC.foldchange(cc) = str2double(t1{4}); % not normalized
    FC.foldchange(cc) = str2double(t1{5}); % normalized
end
fns = fieldnames(FC);
for ii = 1:length(fns)
    FC.(fns{ii}) = FC.(fns{ii})(1:cc);
end


%%   2. Align b/w replicates
% This is done by Alignment.m
% Variable pfit(ci,rr,:)
%   ci: condition
%   rr: replicate
% saved in F:\Greg\PCP_SILAC\Runs\Craig\Craig20170829_tests\Craig20170829_4reps_final_reproduce_filter1/Output/tmp/alignmentfit.mat
fn_pfit = [maindir '/Output/tmp/alignmentfit.mat'];
try
  load(fn_pfit)
  failFlag = 0;
catch
  failFlag = 1;
end

FC.alignedfraction = nan(size(FC.fraction));
for ii = 1:length(FC.protein)
  rep = FC.replicate(ii);
  
  % Hack hack hack hack hack
  % Alignment.m aligns M/L and H/L separately.
  % I calculated FC pre-alignment, though.
  % M/L and H/L alignment should be the same, though.
  % Therefore I'm just averaging them.
  if failFlag==0
    b = nanmean(squeeze(pfit(:,rep,1))); % intercept
    m = nanmean(squeeze(pfit(:,rep,2))); % slope
  else
    % couldn't find the alignment fit!
    b = 0;
    m = 1;
  end
  % Hack hack hack hack hack
  
  FC.alignedfraction(ii) = FC.fraction(ii)*m + b;
end



%%   3. Combine fold changes across replicates

% make structure pFC, which is broken down by protein
clear pFC
pFC.protein = unique(FC.protein);
maxCount = max(cell2mat(cellfun(@(x) sum(ismember(FC.protein,x)),pFC.protein,'un',0)));
pFC.replicate = nan(length(pFC.protein),maxCount);
pFC.fraction = nan(length(pFC.protein),maxCount);
pFC.alignedfraction = nan(length(pFC.protein),maxCount);
pFC.foldchange = nan(length(pFC.protein),maxCount);
for ii = 1:length(pFC.protein)
    I = find(ismember(FC.protein,pFC.protein{ii}));
    for jj = 1:length(I)
        pFC.replicate(ii,jj) = FC.replicate(I(jj));
        pFC.fraction(ii,jj) = FC.fraction(I(jj));
        pFC.alignedfraction(ii,jj) = FC.alignedfraction(I(jj));
        pFC.foldchange(ii,jj) = FC.foldchange(I(jj));
    end
end

% In the same chromatogram, merge all fold changes within 2.5 fractions of each other.
% Use a greedy sliding window.
windowRange = -2:65;
fnames = fieldnames(pFC);
for ii = 1:length(pFC.protein)
    for jj = 1:4
        I = pFC.replicate(ii,:) == jj;
        Ninwindow = ones(size(windowRange));
        cc = 0;
        while sum(Ninwindow)>0
            cc = cc+1;
            fracs = pFC.fraction(ii,I);
            for kk = 1:length(windowRange)
                center = windowRange(kk);
                I_inwindow = fracs>=center-3.5 & fracs<=center+3.5;
                Ninwindow(kk) = sum(I_inwindow);
            end
            Ninwindow(Ninwindow==1) = 0;
            if sum(Ninwindow)==0; continue; end
            
            [x1,I1] = max(Ninwindow);
            center = windowRange(I1);
            I_inwindow = fracs>=center-3.5 & fracs<=center+3.5;
            % merge these fractions
            fracs2merge = fracs(I_inwindow);
            I2merge = find(ismember(pFC.fraction(ii,:), fracs2merge));
            for kk = 1:length(fnames)
                if strcmp(fnames{kk},'protein'); continue; end
                avg = nanmean(pFC.(fnames{kk})(ii,I2merge));
                pFC.(fnames{kk})(ii,I2merge) = nan;
                pFC.(fnames{kk})(ii,I2merge(1)) = avg;
            end
        end
    end
end


%%   4. Group fold changes for a single protein

% Again, use a greedy sliding window
pFC.group = nan(size(pFC.foldchange));
windowRange = -2:65;
fnames = fieldnames(pFC);
for ii = 1:length(pFC.protein)
    Ninwindow = ones(size(windowRange));
    groupN = 0;
    while sum(isnan(pFC.group(ii,:))) > sum(isnan(pFC.fraction(ii,:)))
        groupN = groupN+1;
        fracs = pFC.alignedfraction(ii,:);
        group = pFC.group(ii,:);
        for kk = 1:length(windowRange)
            center = windowRange(kk);
            I_inwindow = fracs>=center-3.5 & fracs<=center+3.5 & isnan(group);
            Ninwindow(kk) = sum(I_inwindow);
        end
        if sum(Ninwindow)==0; continue; end
        
        [x1,I1] = max(Ninwindow);
        center = windowRange(I1);
        I_inwindow = fracs>=center-3.5 & fracs<=center+3.5 & isnan(group);
        % merge these fractions
        fracs2merge = fracs(I_inwindow);
        I2merge = find(ismember(pFC.alignedfraction(ii,:), fracs2merge));
        pFC.group(ii,I2merge) = groupN;
    end
end

% Each group should only appear once per replicate
% Merge groups that occur twice in a replicate
for ii = 1:length(pFC.protein)
    groups = unique(pFC.group(ii,:));
    groups(isnan(groups)) = [];
    for jj = 1:length(groups)
        Igroup = pFC.group(ii,:) == groups(jj);
        for kk = 1:4
            Irep = pFC.replicate(ii,:) == kk;
            I2merge = find(Igroup&Irep);
            if length(I2merge)>1
                % merge
                pFC.fraction(ii,I2merge(1)) = nanmean(pFC.fraction(ii,I2merge),2);
                pFC.alignedfraction(ii,I2merge(1)) = nanmean(pFC.alignedfraction(ii,I2merge),2);
                pFC.foldchange(ii,I2merge(1)) = nanmean(pFC.foldchange(ii,I2merge),2);
                % remove the rest
                pFC.replicate(ii,I2merge(2:end)) = nan;
                pFC.fraction(ii,I2merge(2:end)) = nan;
                pFC.alignedfraction(ii,I2merge(2:end)) = nan;
                pFC.foldchange(ii,I2merge(2:end)) = nan;
                pFC.group(ii,I2merge(2:end)) = nan;
            end
        end
    end
end

% Re-merge groups that are very close
% In the same chromatogram, merge all fold changes within 2.5 fractions of each other.
% Use a greedy sliding window.
windowRange = -2:65;
fnames = fieldnames(pFC);
for ii = 1:length(pFC.protein)
    for jj = 1:4
        I = pFC.replicate(ii,:) == jj;
        Ninwindow = ones(size(windowRange));
        cc = 0;
        while sum(Ninwindow)>0
            cc = cc+1;
            fracs = pFC.fraction(ii,I);
            for kk = 1:length(windowRange)
                center = windowRange(kk);
                I_inwindow = fracs>=center-3.5 & fracs<=center+3.5;
                Ninwindow(kk) = sum(I_inwindow);
            end
            Ninwindow(Ninwindow==1) = 0;
            if sum(Ninwindow)==0; continue; end
            
            [x1,I1] = max(Ninwindow);
            center = windowRange(I1);
            I_inwindow = fracs>=center-3.5 & fracs<=center+3.5;
            % merge these fractions
            fracs2merge = fracs(I_inwindow);
            I2merge = find(ismember(pFC.fraction(ii,:), fracs2merge));
            pFC.group(ii,I2merge(2:end)) = nan;
            for kk = 1:length(fnames)
                if strcmp(fnames{kk},'protein'); continue; end
                if strcmp(fnames{kk},'group'); continue; end
                avg = nanmean(pFC.(fnames{kk})(ii,I2merge));
                pFC.(fnames{kk})(ii,I2merge) = nan;
                pFC.(fnames{kk})(ii,I2merge(1)) = avg;
            end
        end
    end
end

%% Attempt to fill in "missing" fold changes

% Recalculate all fold changes using current fractions
reps = 1:4;
xfrac = 1:60;
for ii = 1:length(pFC.protein)
    prot = pFC.protein{ii,1};
    groups = unique(pFC.group(ii,:));
    groups(isnan(groups)) = [];
    Iraw = find(ismember(rawdata(1).Proteins, prot)) - 1;
    for jj = 1:length(groups)
        I = pFC.group(ii,:) == groups(jj);        
        avgAlignedFraction = nanmean(pFC.alignedfraction(ii,I));
        for kk = 1:length(reps)
            % is there enough data? (>2 data points)
            rep = reps(kk);
            b = nanmean(squeeze(pfit(:,rep,1))); % intercept
            m = nanmean(squeeze(pfit(:,rep,2))); % intercept
            avgFraction = (avgAlignedFraction - b)/m;
            Ifrac = find(xfrac>avgFraction-2.5 & xfrac<avgFraction+2.5);
            thisIraw = Iraw(rep);
            data_num = rawdata(1).data(thisIraw,Ifrac+1);
            data_den = rawdata(2).data(thisIraw,Ifrac+1);
            fc = data_num ./ data_den;
            if sum(not(isnan(fc)))>1
                % found a fold change!
                % find where to put it
                ia = find(pFC.group(ii,:)==groups(jj) & pFC.replicate(ii,:)==rep);
                if isempty(ia)
                    ia = find(isnan(pFC.fraction(ii,:)) & isnan(pFC.replicate(ii,:)),1,'first');
                end
                pFC.replicate(ii,ia) = rep;
                pFC.fraction(ii,ia) = avgFraction;
                pFC.foldchange(ii,ia) = log2(nanmean(fc));
                pFC.group(ii,ia) = groups(jj);
                pFC.alignedfraction(ii,ia) = avgAlignedFraction;
            end
        end
    end
end


%%   5. Choose which fold change groups to keep

% I have up to 7 fold change groups! Some are probably garbage
% How do I choose the good ones to keep?

% Idea 1: Weight each group by the raw data it represents
xfrac = 1:60;
pFC_filtered = pFC;
for ii = 1:length(pFC.protein)
    prot = pFC.protein{ii};
    Iraw = find(ismember(rawdata(1).Proteins, prot)) - 1;
    groups = unique(pFC.group(ii,:));
    groups(isnan(groups)) = [];
    
    if length(groups)<=2
        % If there are 1 or 2 groups, don't touch em
        continue
    else
        % first pass, weigh all groups
        groupWeight = zeros(size(groups));
        I = find(not(isnan(pFC.group(ii,:))));
        for jj = 1:length(I)
            thisGroup = pFC.group(ii,I(jj));
            rep = pFC.replicate(ii,I(jj));
            IthisGroup = groups==thisGroup;
            frac = pFC.fraction(ii,I(jj));
            Ifrac = xfrac>frac-2.5 & xfrac<frac+2.5;
            w1 = nansum(rawdata(1).data(Iraw(rep),Ifrac));
            w2 = nansum(rawdata(2).data(Iraw(rep),Ifrac));
            groupWeight(IthisGroup) = groupWeight(IthisGroup) + w1 + w2;
        end
        if sum(isnan(groupWeight))>0
            error('yoo!')
        end
        
        % second pass
        % skip for now
        
        % choose the two heaviest groups
        [x1,I] = sort(groupWeight,'descend');
        goodGroups = groups(I(1:2));
        badGroups = groups(not(ismember(groups,goodGroups)));
        
        % remove the bad groups
        Ibad = ismember(pFC.group(ii,:), badGroups);
        pFC_filtered.replicate(ii,Ibad) = nan;
        pFC_filtered.fraction(ii,Ibad) = nan;
        pFC_filtered.alignedfraction(ii,Ibad) = nan;
        pFC_filtered.foldchange(ii,Ibad) = nan;
        pFC_filtered.group(ii,Ibad) = nan;
    end
end

% Normalize fold changes
pFC_filtered.foldchange = pFC_filtered.foldchange - nanmedian(pFC_filtered.foldchange(:));


%%   6. Calculate stats

clear dataOut
dataOut.protein = pFC_filtered.protein;
dataOut.avgfraction = nan(length(pFC_filtered.protein), 2);
dataOut.pvalue = nan(length(pFC_filtered.protein), 2);
dataOut.avgfoldchange = nan(length(pFC_filtered.protein), 2);
dataOut.foldchange = nan(length(pFC_filtered.protein), 8);
for ii = 1:length(pFC_filtered.protein)
    groups = unique(pFC_filtered.group(ii,:));
    groups(isnan(groups)) = [];
    if length(groups)>2
        error('uhhh')
    end
    for jj = 1:length(groups)
        Igroup = pFC_filtered.group(ii,:) == groups(jj);
        % fill the fold changes (data) in the order of the replicates
        data = nan(1,4);
        frac = nan(1,4);
        for kk = 1:4
            Irep = pFC_filtered.replicate(ii,:) == kk;
            I = Igroup & Irep;
            if sum(I)>0
                data(kk) = pFC_filtered.foldchange(ii,I);
                frac(kk) = pFC_filtered.alignedfraction(ii,I);
            end
        end
        if length(data)>4
            error('ummm')
        end
        [x1,dataOut.pvalue(ii,jj)] = ttest(data);
        dataOut.foldchange(ii,(jj-1)*4+1:jj*4) = data;
        dataOut.avgfoldchange(ii,jj) = nanmean(data);
        dataOut.avgfraction(ii,jj) = nanmean(frac);
    end
end


%%   6. Connect dataOut with interaction list

header2add = {'Significant fold change? (p<.05 protA)' ... % 1
    'Significant fold change? (p<.05 protB)' ... % 2
    'p-value (protA peak1)' ... % 3
    'p-value (protA peak2)' ... % 4
    'p-value (protB peak1)' ... % 5
    'p-value (protB peak2)' ... % 6
    'Avg log2 fold change (protA peak1)' ... % 7
    'Avg log2 fold change (protA peak2)' ... % 8
    'Avg log2 fold change (protB peak1)' ... % 9
    'Avg log2 fold change (protB peak2)' ... % 10
    'Aligned elution fraction (protA peak1)' ... % 11
    'Aligned elution fraction (protA peak2)' ... % 12
    'Aligned elution fraction (protB peak1)' ... % 13
    'Aligned elution fraction (protB peak2)' ... % 14
    'Log2 fold change (protA peak1 rep1)' ... % 15
    'Log2 fold change (protA peak1 rep2)' ... % 16
    'Log2 fold change (protA peak1 rep3)' ... % 17
    'Log2 fold change (protA peak1 rep4)' ... % 18
    'Log2 fold change (protA peak2 rep1)' ... % 19
    'Log2 fold change (protA peak2 rep2)' ... % 20
    'Log2 fold change (protA peak2 rep3)' ... % 21
    'Log2 fold change (protA peak2 rep4)' ... % 22
	'Log2 fold change (protB peak1 rep1)' ... % 23
    'Log2 fold change (protB peak1 rep2)' ... % 24
    'Log2 fold change (protB peak1 rep3)' ... % 25
    'Log2 fold change (protB peak1 rep4)' ... % 26
    'Log2 fold change (protB peak2 rep1)' ... % 27
    'Log2 fold change (protB peak2 rep2)' ... % 28
    'Log2 fold change (protB peak2 rep3)' ... % 29
    'Log2 fold change (protB peak2 rep4)'}; % 30

dataOut.forPrint = cell(10^5, length(header2add));

% Just read in interaction list stupidly (all strings)
fnin = 'F:\Greg\PCP_SILAC\Runs\Craig\Craig20170829_tests\Craig20170829_4reps_final_reproduce_filter1/Output/Data/Interactions/Final_Interactions_wProtNames_list_50_precision.csv';
fid = fopen(fnin);
datain.head = strsplit(fgetl(fid),',');
I = cellfun('isempty',datain.head);
datain.head(I) = [];
datain.data = cell(10^5,length(datain.head));
cc = 0;
while not(feof(fid))
    cc = cc+1;
    datain.data(cc,:) = strsplit(fgetl(fid),',','collapsedelimiters',0);
    
    protA = datain.data{cc,2};
    protB = datain.data{cc,3};
    Ia = find(ismember(dataOut.protein,protA));
    Ib = find(ismember(dataOut.protein,protB));
    
    % is there a significant change?
    % 1. protA
    I = dataOut.pvalue(Ia,:)<.05;
    if sum(I)>0
        tmp = dataOut.avgfoldchange(Ia,I);
        tmp(tmp==0) = [];
        if sum(tmp>0)==length(tmp)
            % all increase?
            dataOut.forPrint{cc,1} = 'Increase';
        elseif sum(tmp<0)==length(tmp)
            % all decrease
            dataOut.forPrint{cc,1} = 'Decrease';
        else
            % both?
            dataOut.forPrint{cc,1} = 'Increase/Decrease';
        end
    end
    
    % 2. protB
    I = dataOut.pvalue(Ib,:)<.05;
    if sum(I)>0
        tmp = dataOut.avgfoldchange(Ib,I);
        tmp(tmp==0) = [];
        if sum(tmp>0)==length(tmp)
            % all increase?
            dataOut.forPrint{cc,2} = 'Increase';
        elseif sum(tmp<0)==length(tmp)
            % all decrease
            dataOut.forPrint{cc,2} = 'Decrease';
        else
            % both?
            dataOut.forPrint{cc,2} = 'Increase/Decrease';
        end
    end
    
    % p-values
    dataOut.forPrint{cc,3} = dataOut.pvalue(Ia,1);
    dataOut.forPrint{cc,4} = dataOut.pvalue(Ia,2);
    dataOut.forPrint{cc,5} = dataOut.pvalue(Ib,1);
    dataOut.forPrint{cc,6} = dataOut.pvalue(Ib,2);
    
    % avg fold change
    dataOut.forPrint{cc,7} = dataOut.avgfoldchange(Ia,1);
    dataOut.forPrint{cc,8} = dataOut.avgfoldchange(Ia,2);
    dataOut.forPrint{cc,9} = dataOut.avgfoldchange(Ib,1);
    dataOut.forPrint{cc,10} = dataOut.avgfoldchange(Ib,2);
    
    % Aligned elution peaks
    dataOut.forPrint{cc,11} = dataOut.avgfraction(Ia,1);
    dataOut.forPrint{cc,12} = dataOut.avgfraction(Ia,2);
    dataOut.forPrint{cc,13} = dataOut.avgfraction(Ib,1);
    dataOut.forPrint{cc,14} = dataOut.avgfraction(Ib,2);
    
    % Fold changes
    ccount = 14;
    for kk = 1:2
        for jj = 1:4
            ccount = ccount+1;
            dataOut.forPrint{cc,ccount} = dataOut.foldchange(Ia,jj + (kk-1)*4);
        end
    end
    for kk = 1:2
        for jj = 1:4
            ccount = ccount+1;
            dataOut.forPrint{cc,ccount} = dataOut.foldchange(Ib,jj + (kk-1)*4);
        end
    end
end
datain.data = datain.data(1:cc,:);
dataOut.forPrint = dataOut.forPrint(1:cc,:);
fclose all;


%%   6. Write fold changes to final_interactions.csv

fnout = [maindir '/Output/Data/Interactions/Final_Interactions_wProtNames_FoldChanges_list_' ...
  numstr(round(user.desiredPrecision*100)) '_precision.csv'];
fid = fopen(fnout,'w');
% header
for ii = 1:length(datain.head)
    fprintf(fid,'%s,', datain.head{ii});
end
for ii = 1:length(header2add)
    fprintf(fid,'%s,', header2add{ii});
end
fprintf(fid,'\n');
for ii = 1:cc
    % print interaction list data
    for jj = 1:size(datain.data,2)
        fprintf(fid,'%s,',datain.data{ii,jj});
    end
    
    % protein A
    for jj = 1:size(dataOut.forPrint,2)
        if jj>2
            fprintf(fid,'%6.4f,', dataOut.forPrint{ii,jj});
        else
            fprintf(fid,'%s,', dataOut.forPrint{ii,jj});
        end
    end
    fprintf(fid,'\n');
end
fclose all;



