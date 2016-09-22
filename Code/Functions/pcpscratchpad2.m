% May 31
%
% Think I just fixed THE BUG, and that it had to do with my Corum files + TP_Matrix error.

%% Once I reduce to just A-B interactions, are my corum files the same as Nick's?

% Apoptosis data
ff{1} = '/Users/Mercy/Academics/Foster/NickCodeData/Old runs/GregPCP_20160517/Data/Corum_correctly_formated_Uniprot_IDs.csv';
ff{2} = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/Output/tmp/Corum_pairwise.csv';
corList_apt = cell(2,1);
for ii = 1:2
  fid=fopen(ff{ii}, 'rt');    %Input corum data base.
  Corum_Import= textscan(fid,'%s\t', 'Delimiter',',');
  fclose(fid);
  No=length(Corum_Import{1})/2;
  tmp = reshape(Corum_Import{1,1},2,No)';
  corList_apt{ii} = cell(size(tmp,1),1);
  for jj = 1:size(tmp,1)
    tmp(jj,:) = sort(tmp(jj,:));
    corList_apt{ii}{jj} = strjoin(tmp(jj,:),'_');
  end
  
  % reduce to unique
  %[~,idx]=unique(cell2mat(corList_apt{ii}),'rows');
  %corList_apt{ii} =  corList_apt{ii}(idx,:);
  corList_apt{ii} = unique(corList_apt{ii});
end
I = intersect(corList_apt{1},corList_apt{2});
figure,hold on
myVenn2([length(corList_apt{1}) length(corList_apt{2})],length(I))

%%
% Tissue data
ff{1} = '/Users/Mercy/Academics/Foster/Tissue_PCPSILAC/PCPSILAC_analysis/Input/Mapped_mouse_Corum_list_20150109.csv';
ff{2} = '/Users/Mercy/Academics/Foster/Tissue_PCPSILAC/PCPSILAC_analysis/Output/tmp/Corum_pairwise.csv';
corList_tis = cell(2,1);
for ii = 1:2
  fid=fopen(ff{ii}, 'rt');    %Input corum data base.
  Corum_Import= textscan(fid,'%s\t', 'Delimiter',',');
  fclose(fid);
  No=length(Corum_Import{1})/2;
  tmp = reshape(Corum_Import{1,1},2,No)';
  corList_tis{ii} = cell(size(tmp,1),1);
  for jj = 1:size(tmp,1)
    tmp(jj,:) = sort(tmp(jj,:));
    corList_tis{ii}{jj} = strjoin(tmp(jj,:),'_');
  end
  
  % reduce to unique
  %[~,idx]=unique(cell2mat(corList_apt{ii}),'rows');
  %corList_apt{ii} =  corList_apt{ii}(idx,:);
  corList_tis{ii} = unique(corList_tis{ii});
end
I = intersect(corList_tis{1},corList_tis{2});
figure,hold on
myVenn2([length(corList_tis{1}) length(corList_tis{2})],length(I))



%% Classify on all replicates together
% use Apoptosis data
% run pcpsilac_apoptosis first

Proteins = cell(6,1);
Dists = cell(6,1);
uniqueProteins = [];
uniqueProteins_MajorIDs = [];
for ii = 3:3
  sf = [user.maindir '/Output/tmp/' 'score_rep' num2str(ii) '.mat'];
  load(sf)
  Proteins{ii} = Protein;
  Dists{ii} = Dist;
  clear scoreMatrix TP_Matrix possibleInts Protein inverse_self Chromatograms Dist
  
  [uniqueProteins,I1,I2] = unique([uniqueProteins; Proteins{ii}.NoIsoform]);
  for jj = 1:length(uniqueProteins)
    I = find(ismember(Proteins{ii}.NoIsoform,uniqueProteins{jj}));
    I2 = find(~cellfun('isempty',Proteins{ii}.MajorID_NoIsoforms(I,:)));
    for kk = 1:length(I2)
      uniqueProteins_MajorIDs{jj,I2(kk)} = Proteins{ii}.MajorID_NoIsoforms{I,I2(kk)};
    end
  end
  %   I1 = ismember(Proteins{ii}.Isoform, uniqueProteins); %where are these proteins in Protein?
  %   I2 = ismember(uniqueProteins, Proteins{ii}.Isoform); %where are these proteins in uniqueProteins?
  %   uniqueProteins_MajorIDs(I2,:) = Proteins{ii}.MajorID_NoIsoforms(I1,:);
end


% Make Dist by merging replicates somehow
clear Dist
for ii = 3:3
  I = find(ismember(uniqueProteins,Proteins{ii}.NoIsoform));
  [~,I2] = unique(Proteins{ii}.NoIsoform);
  fn = fieldnames(Dists{ii});
  for jj = 1:length(fn)
    fn1 = [fn{jj} '_' num2str(ii)];
    Dist.(fn1) = nan(length(uniqueProteins),length(uniqueProteins));
    Dist.(fn1)(I,I) = Dists{ii}.(fn{jj})(I2,I2);
    
    % remove the field to save space
    Dists{ii}.(fn{jj}) = [];
  end
end
clear Dists


% Make TP_Matrix
% Corum binary interactions
fid=fopen(user.corumpairwisefile, 'rt');    %Input corum data base.
Corum_Import= textscan(fid,'%s\t', 'Delimiter',',');
fclose(fid);
No=length(Corum_Import{1})/2;
Corum_Protein_names = (reshape(Corum_Import{1,1},2,No)');
Unique_Corum = unique(Corum_Protein_names);
Pos_Corum_proteins = Unique_Corum(ismember(Unique_Corum, uniqueProteins));
% Make Corum_Dataset, which is all interactions that are in corum and our dataset
cc = 0;
Corum_Dataset = cell(1000,1);
for jj = 1:length(Corum_Protein_names); %Write out values as i
  Prot1 = Corum_Protein_names(jj,1);
  Prot2 = Corum_Protein_names(jj,2);
  
  % find these protein names in our sample
  tmp1 = sum(strcmp(Prot1,Pos_Corum_proteins));
  tmp2 = sum(strcmp(Prot2,Pos_Corum_proteins));
  
  if tmp1>0 && tmp2>0
    cc = cc+1;
    Corum_Dataset{cc} = [Prot1{1},'-',Prot2{1}];
  end
end
TP_Matrix= zeros(length(uniqueProteins),length(uniqueProteins));
%Corum_Dataset_expanded = cell(size(Corum_Dataset));
ii = 0;
for cc=1:size(Corum_Dataset,1)
  I_hyphen = strfind(Corum_Dataset{cc},'-');
  Prot1 = Corum_Dataset{cc}(1:I_hyphen-1);
  Prot2 = Corum_Dataset{cc}(I_hyphen+1:end);
  
  % find where Prot1 and Prot2 are in Protein.MajorID_NoIsoforms
  [i1, ~] = find(strcmp(Prot1,uniqueProteins_MajorIDs));
  [i2, ~] = find(strcmp(Prot2,uniqueProteins_MajorIDs));
  
  TP_Matrix(i1,i2) = 1;
  TP_Matrix(i2,i1) = 1;
end

% Make possibleInts
% check which individual proteins are in corum
% Important! Also check whether any member of that protein's group is in corum
Ngroup = size(uniqueProteins_MajorIDs,2);
inCorum = zeros(Ngroup,length(uniqueProteins));
for ni = 1:Ngroup
  tmp = uniqueProteins_MajorIDs(:,ni);
  % ismember can't handle empty cells, so as a workaround fill with a string that won't be in Corum
  isemptyCellArray = cellfun('isempty',tmp);
  tmp(isemptyCellArray) = {'--fake'};
  inCorum(ni,:) = ismember(tmp,Unique_Corum);
end
inCorum = sum(inCorum,1);
% Calculate possible interaction matrix, asks are both proteins in corum
Int_matrix = inCorum'*inCorum;
% Create self interaction matrix
Self_Int_matrix = zeros(size(TP_Matrix));
for ii = 1:length(uniqueProteins)
  for jj = 1:length(uniqueProteins)
    if isequal(uniqueProteins{ii},uniqueProteins{jj})
      Self_Int_matrix(ii,jj)=1;
    end
  end
end
inverse_self = ~Self_Int_matrix; %creation of inverse of self interactions
clear Self_Int_matrix
% Possible interactions, i.e. both in corum and not a self interaction
possibleInts = inverse_self & Int_matrix;


% Classify!
[scoreMatrix,feats] = scorenb(Dist,possibleInts,TP_Matrix);
score = nanmedian(scoreMatrix,2);


% Make PR curve
x = [0:0.1:0.5 0.55:0.05:0.9 0.91:.01:.99 .991:.001:.999 .9991:.0001:.9999 .99991:.00001:.99999 .999991:.000001:.999999...
  .9999991:.0000001:.9999999 .99999991:.00000001:.99999999];
Rec = zeros(size(x));
Prec = zeros(size(x));
Nint = zeros(size(x));
for xi = 1:length(x)
  xi
  TP = sum(score(:)>x(xi) & TP_Matrix(:)==1 & possibleInts(:)==1);
  FP = sum(score(:)>x(xi) & TP_Matrix(:)==0 & possibleInts(:)==1);
  FN = sum(score(:)<x(xi) & TP_Matrix(:)==1 & possibleInts(:)==1);
  
  Prec(xi) = TP / (TP + FP);
  Rec(xi) = TP / (TP + FN);
  Nint(xi) = sum(score(:)>x(xi));
end
figure,plot(Nint,Prec)



%% Classify on all replicates together, v0.2
% Don't make TP_Matrix, possibleInts.
% Cobble them together from the save files


% make uniqueProteins
disp('Making uniqueProteins...')
Proteins = cell(6,1);
uniqueProteins = [];
uniqueProteins_MajorIDs = [];
for ii = 1:6
  sf = [user.maindir '/Output/tmp/' 'score_rep' num2str(ii) '.mat'];
  load(sf)
  uniqueProteins = unique([uniqueProteins; Protein.Isoform]);
  Proteins{ii} = Protein;
  clear scoreMatrix TP_Matrix possibleInts Protein inverse_self Chromatograms Dist
end


disp('Making TP_Matrix and possibleInts...')
TP_Matrix2 = zeros(length(uniqueProteins),length(uniqueProteins));
possibleInts2 = zeros(length(uniqueProteins),length(uniqueProteins));
clear Dist2
for ii = 1:6
  ii
  sf = [user.maindir '/Output/tmp/' 'score_rep' num2str(ii) '.mat'];
  load(sf)
  
  Iunq = zeros(size(Proteins{ii}.Isoform));
  for jj = 1:length(Protein.Isoform)
    Iunq(jj) = find(ismember(uniqueProteins,Proteins{ii}.Isoform{jj}));
  end
  
  for jj = 1:length(Protein.Isoform)
    for kk = 1:length(Protein.Isoform)
      TP_Matrix2(Iunq(jj),Iunq(kk)) = TP_Matrix(jj,kk);
      TP_Matrix2(Iunq(kk),Iunq(jj)) = TP_Matrix(kk,jj);
      
      possibleInts2(Iunq(jj),Iunq(kk)) = possibleInts(jj,kk);
      possibleInts2(Iunq(kk),Iunq(jj)) = possibleInts(kk,jj);
    end
  end
  
  fn = fieldnames(Dist);
  for mm = 1:length(fn)
    mm
    fn1 = [fn{mm} '_' num2str(ii)];
    A = nan(length(uniqueProteins),length(uniqueProteins));
    B = Dist.(fn{mm});
    for jj = 1:length(Protein.Isoform)
      for kk = 1:length(Protein.Isoform)
        A(Iunq(jj),Iunq(kk)) = B(jj,kk);
        A(Iunq(kk),Iunq(jj)) = B(kk,jj);
      end
    end
    Dist2.(fn1) = A;
  end
  clear scoreMatrix TP_Matrix possibleInts Protein inverse_self Chromatograms Dist A B
end
Dist2

% Classify!
disp('Classifying...')
[scoreMatrix,feats] = scorenb(Dist2,possibleInts2,TP_Matrix2);
score = nanmedian(scoreMatrix,2);
class = nan(size(score(:)));
class(possibleInts2(:)==1) = TP_Matrix2(possibleInts2(:)==1);


% Make PR curve
disp('Making PR curve...')
[~, vv] = calcPPIthreshold(score, class, .8);
vv.Ninteract = zeros(size(vv.scoreRange));
for ii = 1:length(vv.scoreRange)
  vv.Ninteract(ii) = sum(score > vv.scoreRange(ii));
end

mySound

figure
subplot(2,1,1)
semilogx(vv.Ninteract,vv.precRange),grid on
title(user.maindir)
ylabel('Precision')
xlabel('Number of interactions')
subplot(2,1,2)
plot(vv.recRange,vv.precRange),grid on
ylabel('Precision')
xlabel('Recall')



%% makeFigures_comparison is broken for Apoptosis + individual protein figures
% run makeFigures_comparison to the first figure

ii = 1;
rep = 3;

% get replicate + protein name
protName = Finalised_Master_Gaussian_list.Protein_name{ii,1};
%rep_protName = Finalised_Master_Gaussian_list.Replicate_Protein_identifier{ii};
rep_protName = [num2str(rep) '_' protName]

% Find this replicate + protein in Combined_Gaussians
Icg = find(Combined_Gaussians.Replicate == rep & strcmp(protName,Combined_Gaussians.Protein_name));
if isempty(Icg)
  continue
end
[C_comp,Iccomp] = sort(Combined_Gaussians.Center(Icg));

% Get channel names and indices
chan_den = Combined_Gaussians.denominatorChannnel{Icg};
chan_num = Combined_Gaussians.numeratorChannnel{Icg};
Iden = find(ismember(user.silacratios,chan_den));
Inum = find(ismember(user.silacratios,chan_num));

% get raw data
Iraw = find(ismember(txt_val{Iden}(:,1), rep_protName));
raw_num = num_val{Inum}(Iraw,2:end);
raw_den = num_val{Iden}(Iraw,2:end);

% get gaussian params
Igd_den = find(ismember(GaussData{Iden}(:,2),rep_protName));
hd = (GaussData{Iden}(Igd_den,7));
cd = (GaussData{Iden}(Igd_den,8));
wd = (GaussData{Iden}(Igd_den,9));
gp_den = zeros(length(hd),3);
for jj = 1:length(hd)
  gp_den(jj,1) = str2num(hd{jj});
  gp_den(jj,2) = str2num(cd{jj});
  gp_den(jj,3) = str2num(wd{jj});
end
Igd_num = find(ismember(GaussData{Inum}(:,2),rep_protName));
hn = (GaussData{Inum}(Igd_num,7));
cn = (GaussData{Inum}(Igd_num,8));
wn = (GaussData{Inum}(Igd_num,9));
gp_num = zeros(length(hn),3);
for jj = 1:length(hn)
  gp_num(jj,1) = str2num(hn{jj});
  gp_num(jj,2) = str2num(cn{jj});
  gp_num(jj,3) = str2num(wn{jj});
end

% make gaussian curves
xfit = linspace(0,length(raw_num)+1,101);
yfit_den = zeros(size(xfit));
for jj = 1:length(cd)
  yfit_den = yfit_den + gp_den(jj,1)*exp(-((xfit-gp_den(jj,2))/gp_den(jj,3)).^2);
end
xfit = linspace(0,length(raw_num)+1,101);
yfit_num = zeros(size(xfit));
for jj = 1:length(cn)
  yfit_num = yfit_num + gp_num(jj,1)*exp(-((xfit-gp_num(jj,2))/gp_num(jj,3)).^2);
end

% get fold changes for each gaussian
fold_raw = Combined_Gaussians.log2_of_gaussians(Icg);
fold_norm = Combined_Gaussians.log2_normalised_gaussians(Icg);

figure
subplot(3,1,1), hold on % raw chromatograms
xraw = 1:length(raw_num);
scatter(xraw,raw_num,20,colour_to_use(2,:),'filled')
scatter(xraw,raw_den,20,colour_to_use(1,:),'filled')
plot(xraw,raw_num,'color',colour_to_use(2,:));
plot(xraw,raw_den,'color',colour_to_use(1,:));
ax = axis;

subplot(3,1,2), hold on % fit gaussians, numerator channel
plot(xfit,yfit_num,'color',colour_to_use(2,:));
plot(xfit,yfit_den,'color',colour_to_use(1,:));
axis(ax)



%% try it with params from Finalised_Master_Gaussian_list

ii = 2;
rep = 3;

% get replicate + protein name
protName = Finalised_Master_Gaussian_list.Protein_name{ii,1};
rep_protName = [num2str(rep) '_' protName];

% get gaussian params from Finalised_Master_Gaussian_list
I1 = Finalised_Master_Gaussian_list.Replicate(ii,:) == rep & strcmp(Finalised_Master_Gaussian_list.Channel(ii,:), user.silacratios(1));
I2 = Finalised_Master_Gaussian_list.Replicate(ii,:) == rep & strcmp(Finalised_Master_Gaussian_list.Channel(ii,:), user.silacratios(2));
C1 = Finalised_Master_Gaussian_list.Center(ii,I1);
H1 = Finalised_Master_Gaussian_list.Height(ii,I1);
W1 = Finalised_Master_Gaussian_list.Width(ii,I1);
C2 = Finalised_Master_Gaussian_list.Center(ii,I2);
H2 = Finalised_Master_Gaussian_list.Height(ii,I2);
W2 = Finalised_Master_Gaussian_list.Width(ii,I2);

% get raw data
Iraw = find(ismember(txt_val{1}(:,1), rep_protName));
raw1 = num_val{1}(Iraw,2:end);
raw2 = num_val{2}(Iraw,2:end);

% make gaussian curves
xfit = linspace(0,length(raw1)+1,101);
yfit1 = zeros(size(xfit));
yfit2 = zeros(size(xfit));
for ii = 1:length(C1)
  yfit1 = yfit1 + H1(ii) * exp(-(xfit - C1(ii)).^2 / W1(ii)^2);
end
for ii = 1:length(C2)
  yfit2 = yfit2 + H2(ii) * exp(-(xfit - C2(ii)).^2 / W2(ii)^2);
end

figure,hold on
plot(raw1,'r')
plot(raw2,'g')
legend([user.silacratios{1} '1'],[user.silacratios{2} '2'],'location','northwest')
ax = axis;
plot(xfit,yfit1,'--r')
plot(xfit,yfit2,'--g')



%% try it with params from GaussData

ii = 4;
rep = 3;

% get replicate + protein name
protName = Finalised_Master_Gaussian_list.Protein_name{ii,1};
rep_protName = [num2str(rep) '_' protName]

% get gaussian params from GaussData
I1 = find(strcmp(GaussData{2}(:,2), rep_protName));
I2 = find(strcmp(GaussData{1}(:,2), rep_protName));
clear H1 C1 W1 H2 C2 W2
for ii = 1:length(I1)
  H1(ii) = str2num(GaussData{2}{I1(ii),7});
  C1(ii) = str2num(GaussData{2}{I1(ii),8});
  W1(ii) = str2num(GaussData{2}{I1(ii),9});
end
% for ii = 1:length(I2)
%   H2(ii) = str2num(GaussData{1}{I2(ii),7});
%   C2(ii) = str2num(GaussData{1}{I2(ii),8});
%   W2(ii) = str2num(GaussData{1}{I2(ii),9});
% end
H2 = [2.0055 2.8008 2.8012 3.4078];
C2 = [17.5346 30.5027 39.2401 46.7942];
W2 = [2.9793 5.4812 5.3184 2.5040];

% get raw data
Iraw = find(ismember(txt_val{1}(:,1), rep_protName));
raw1 = num_val{1}(Iraw,2:end);
Iraw = find(ismember(txt_val{2}(:,1), rep_protName));
raw2 = num_val{2}(Iraw,2:end);

% make gaussian curves
xfit = linspace(0,length(raw1)+1,101);
yfit1 = zeros(size(xfit));
yfit2 = zeros(size(xfit));
for ii = 1:length(C1)
  yfit1 = yfit1 + H1(ii) * exp(-(xfit - C1(ii)).^2 / W1(ii)^2);
end
for ii = 1:length(C2)
  yfit2 = yfit2 + H2(ii) * exp(-(xfit - C2(ii)).^2 / W2(ii)^2);
end

figure,hold on
plot(raw1,'r')
plot(raw2,'g')
legend([user.silacratios{1} '1'],[user.silacratios{2} '2'],'location','northwest')
ax = axis;
plot(xfit,yfit1,'--r')
plot(xfit,yfit2,'--g')



%%

sf = '/Users/Mercy/Downloads/Final_Interactions_list_50_precision.csv';
fid = fopen(sf);
head = fgetl(fid);
data = zeros(10^5,5);
cc = 0;
while ~feof(fid)
  t = fgetl(fid);
  t1 = strsplit(t,',');
  cc = cc+1;
  data(cc,1) = str2num(t1{11});
  data(cc,2) = str2num(t1{12});
  data(cc,3) = str2num(t1{13});
  data(cc,4) = str2num(t1{16});
  data(cc,5) = str2num(t1{17});
  if mod(cc,1000)==0;disp(num2str(cc));end
end
fclose(fid);



%% Try out new Comparison method

% 0. Run up to Load Data in Gauss_Build for apoptosis

% 1. Choose protein, plot it
%for ii = 2:size(txt_val{1},1)
ii = 4;
protName = txt_val{1}{ii,1};

I1 = strfind(txt_val{1}(:,1),protName);
I1 = find(~cellfun('isempty', I1)) - 1;
I2 = strfind(txt_val{2}(:,1),protName);
I2 = find(~cellfun('isempty', I2)) - 1;

y1 = rawdata{1}(I1,:);
y2 = rawdata{2}(I2,:);
x1 = repmat(1:55,3,1);
x2 = repmat(1:55,3,1);
I = ~isnan(y1(:));
y1 = y1(I);
x1 = x1(I);
I = ~isnan(y2(:));
y2 = y2(I);
x2 = x2(I);

%pause
%close all
%end


% 2. Try a t-test(log(data))
xi = 25;
i1 = x1 == xi;
i2 = x1 == xi;
p = ttest3(log(y1(i1)),log(y2(i2)));


% 3. Get CI from single fit

% Define fit type
ft = fittype('gauss1');

% Define fit options
LB = [0.2 0 1]; % H, C, W
UB = [inf 55 inf]; % H, C, W
fo = fitoptions('method','NonlinearLeastSquares','Robust','LAR','Lower',repmat(LB,1,1),...
  'Upper',repmat(UB,1,1),'MaxIter',400,'Display','off');

[ics1(1), I] = max(y1);
ics1(2) = x1(I);
ics1(3) = 2;
ics1 = ics1 + (rand(1,3)-.5) .* [1 3 1];
[ics2(1), I] = max(y2);
ics2(2) = x2(I);
ics2(3) = 2;
ics2 = ics2 + (rand(1,3)-.5) .* [1 3 1];

set(fo,'startpoint',ics1); % Include ICs in fit options
curveFit1 = fit(x1,y1,ft,fo);
CI1 = confint(curveFit1);
set(fo,'startpoint',ics2); % Include ICs in fit options
curveFit2 = fit(x2,y2,ft,fo);
CI2 = confint(curveFit2);

x = 1:55;
y1a = CI1(1,1) * exp(- ((x - CI1(1,2))/CI1(1,3)).^2);
y1b = CI1(2,1) * exp(- ((x - CI1(2,2))/CI1(2,3)).^2);
y2a = CI2(1,1) * exp(- ((x - CI2(1,2))/CI2(1,3)).^2);
y2b = CI2(2,1) * exp(- ((x - CI2(2,2))/CI2(2,3)).^2);

I_reverse = sort(1:length(x),'descend');

figure,hold on
scatter(x1,y1,30,'r','filled')
scatter(x2,y2,30,'b','filled')
legend('Condition 1', 'Condition 2')
patch([x x(I_reverse)], [y1a y1b(I_reverse)],[1 .9 .9],'edgecolor',[1 .9 .9])
patch([x x(I_reverse)], [y2a y2b(I_reverse)],[.9 .9 1],'edgecolor',[.9 .9 1])
scatter(x1,y1,30,'r','filled')
scatter(x2,y2,30,'b','filled')
set(gca,'xtick',[14 28 38])
grid on



% % 3. Get CI from bootstrap
%
% % Define fit type
% ft = fittype('gauss1');
%
% % Define fit options
% LB = [0.2 0 1]; % H, C, W
% UB = [inf 55 inf]; % H, C, W
% fo = fitoptions('method','NonlinearLeastSquares','Robust','LAR','Lower',repmat(LB,1,1),...
%   'Upper',repmat(UB,1,1),'MaxIter',400,'Display','off');
%
% iterMax = 100;
% fit1_b = zeros(iterMax,length(x));
% fit2_b = zeros(iterMax,length(x));
% for iter = 1:iterMax
%   iter
%   i1 = randsample(length(y1),length(y1),1);
%   i2 = randsample(length(y2),length(y2),1);
%
%   y1_b = y1(i1);
%   x1_b = x1(i1);
%   y2_b = y2(i2);
%   x2_b = x2(i2);
%
%   [ics1(1), I] = max(y1_b);
%   ics1(2) = x1_b(I);
%   ics1(3) = 2;
%   ics1 = ics1 + (rand(1,3)-.5) .* [1 3 1];
%   [ics2(1), I] = max(y2_b);
%   ics2(2) = x2_b(I);
%   ics2(3) = 2;
%   ics2 = ics2 + (rand(1,3)-.5) .* [1 3 1];
%
%   set(fo,'startpoint',ics1); % Include ICs in fit options
%   curveFit1 = fit(x1_b,y1_b,ft,fo);
%   cf1 = coeffvalues(curveFit1);
%   set(fo,'startpoint',ics2); % Include ICs in fit options
%   curveFit2 = fit(x2_b,y2_b,ft,fo);
%   cf2 = coeffvalues(curveFit2);
%
%   fit1_b(iter,:) = cf1(1) * exp(- ((x - cf1(2))/cf1(3)).^2);
%   fit2_b(iter,:) = cf2(1) * exp(- ((x - cf2(2))/cf2(3)).^2);
% end
% y1a_b = mean(fit1_b) + std(fit1_b);
% y1b_b = mean(fit1_b) - std(fit1_b);
% y2a_b = mean(fit2_b) + std(fit2_b);
% y2b_b = mean(fit2_b) - std(fit2_b);
%
% I_reverse = sort(1:length(x),'descend');
%
% figure,hold on
% scatter(x1,y1,30,'r','filled')
% scatter(x2,y2,30,'b','filled')
% legend('Condition 1', 'Condition 2')
% patch([x x(I_reverse)], [y1a_b y1b_b(I_reverse)],[1 .9 .9],'edgecolor',[1 .9 .9])
% patch([x x(I_reverse)], [y2a_b y2b_b(I_reverse)],[.9 .9 1],'edgecolor',[.9 .9 1])
% plot(x,mean(y1_b))
% scatter(x1,y1,30,'r','filled')
% scatter(x2,y2,30,'b','filled')
% set(gca,'xtick',[14 28 38])
% grid on
%


% 3. Get CI from resampling
iterMax = 1000;
fit1 = zeros(iterMax,1);
fit2 = zeros(iterMax,1);
for iter = 1:iterMax
  % choose params1 from CI
  cf1 = zeros(1,size(CI1,2));
  cf2 = zeros(1,size(CI1,2));
  for ci = 1:size(CI1,2)
    mu = mean(CI1(:,ci));
    sd = abs(CI1(1,ci) - mu) / 1.96;
    cf1(ci) = randn(1)*sd + mu;
    
    mu = mean(CI2(:,ci));
    sd = abs(CI2(1,ci) - mu) / 1.96;
    cf2(ci) = randn(1)*sd + mu;
  end
  
  fit1(iter) = cf1(1) * exp(- ((xi - cf1(2))/cf1(3)).^2);
  fit2(iter) = cf2(1) * exp(- ((xi - cf2(2))/cf2(3)).^2);
end


%% Try the above (new Comparison method) for a lot of proteins

% 0. Run up to Load Data in Gauss_Build for apoptosis

% 1. Choose protein, plot it
%pttest = nan(100,1);
%pboot = nan(100,1);
%absD = nan(100,1);
for ii = 101:1000%size(txt_val{1},1)
  protName = txt_val{1}{ii,1};
  
  I1 = strfind(txt_val{1}(:,1),protName);
  I1 = find(~cellfun('isempty', I1)) - 1;
  I2 = strfind(txt_val{2}(:,1),protName);
  I2 = find(~cellfun('isempty', I2)) - 1;
  
  y1 = rawdata{1}(I1,:);
  y2 = rawdata{2}(I2,:);
  x1 = repmat(1:55,3,1);
  x2 = repmat(1:55,3,1);
  I = ~isnan(y1(:));
  y1 = y1(I);
  x1 = x1(I);
  I = ~isnan(y2(:));
  y2 = y2(I);
  x2 = x2(I);
  
  if length(x1)<5 || length(x2)<5
    s = [protName ': Abort'];
    disp(s)
    continue;
  end
  
  % 2. Try a t-test(log(data))
  [~,I] = max(y1);
  xi = x1(I);
  xi = xi(1) + [-2 -1 0 1 2];
  xi(xi<1 | xi>55) = [];
  i1 = x1>=min(xi) & x1<=max(xi);
  i2 = x2>=min(xi) & x2<=max(xi);
  data1 = nan(length(xi),3);
  data2 = nan(length(xi),3);
  for jj = 1:length(xi)
    ii1 = x1==xi(jj);
    data1(jj,1:sum(ii1)) = y1(ii1);
    ii2 = x2==xi(jj);
    data2(jj,1:sum(ii2)) = y2(ii2);
  end
  data1 = nanmean(data1);
  data2 = nanmean(data2);
  data1(isnan(data1)) = [];
  data2(isnan(data2)) = [];
  pttest(ii) = ttest3(log(data1),log(data2));
  absD(ii) = (nanmean(y1(i1(:))) - nanmean(y2(i2(:))));
  
  
  % 3. Get CI from single fit
  
  % Define fit type
  ft = fittype('gauss1');
  
  % Define fit options
  LB = [0 -inf 0]; % H, C, W
  UB = [inf inf inf]; % H, C, W
  fo = fitoptions('method','NonlinearLeastSquares','Robust','LAR','Lower',repmat(LB,1,1),...
    'Upper',repmat(UB,1,1),'MaxIter',400,'Display','off');
  
  [ics1(1), I] = max(y1);
  ics1(2) = x1(I);
  ics1(3) = 2;
  ics1 = ics1 + (rand(1,3)-.5) .* [1 3 1];
  [ics2(1), I] = max(y2);
  ics2(2) = x2(I);
  ics2(3) = 2;
  ics2 = ics2 + (rand(1,3)-.5) .* [1 3 1];
  
  set(fo,'startpoint',ics1); % Include ICs in fit options
  curveFit1 = fit(x1,y1,ft,fo);
  CI1 = confint(curveFit1);
  set(fo,'startpoint',ics2); % Include ICs in fit options
  curveFit2 = fit(x2,y2,ft,fo);
  CI2 = confint(curveFit2);
  
  
  
  % 3. Get CI from resampling
  iterMax = 1000;
  fit1 = zeros(iterMax,1);
  fit2 = zeros(iterMax,1);
  for iter = 1:iterMax
    % choose params1 from CI
    cf1 = zeros(1,size(CI1,2));
    cf2 = zeros(1,size(CI1,2));
    for ci = 1:size(CI1,2)
      mu = mean(CI1(:,ci));
      sd = abs(CI1(1,ci) - mu) / 1.96;
      cf1(ci) = randn(1)*sd + mu;
      
      mu = mean(CI2(:,ci));
      sd = abs(CI2(1,ci) - mu) / 1.96;
      cf2(ci) = randn(1)*sd + mu;
    end
    
    fit1(iter) = nanmean(cf1(1) * exp(- ((xi - cf1(2))/cf1(3)).^2));
    fit2(iter) = nanmean(cf2(1) * exp(- ((xi - cf2(2))/cf2(3)).^2));
  end
  
  p1 = nansum(fit1-fit2 < 0) / iterMax;
  p2 = nansum(fit2-fit1 < 0) / iterMax;
  pboot(ii) = min([p1 p2]);
  
  s = [protName ': ' num2str(absD(ii)) ', ' num2str(pttest(ii)) ', ' num2str(pboot(ii))];
  disp(s)
  
end



%% Look for Mike Carlson's gold standard TP

sf = [maindir '/Output/tmp/' 'score_rep' num2str(1) '.mat'];
load(sf)
scoreMatrix = reshape(scoreMatrix,sqrt(length(scoreMatrix)), sqrt(length(scoreMatrix)));

nameGold = {'b0408' 'b0409';'b0177' 'b2512'; 'b0722' 'b0723'; 'b0722' 'b0724'; 'b0723' 'b0724'};

tmp = cell(size(Unique_MajorProteinIDs,1),1);
for ii = 1:length(tmp)
  tmp{ii} = Unique_MajorProteinIDs(ii,:);
end

Igold = zeros(size(nameGold,1),2);
for ii = 1:size(nameGold,1)
  I1 = find(ismember(tmp, nameGold{ii,1}));
  I2 = find(ismember(tmp, nameGold{ii,2}));
  if ~isempty(I1)
    Igold(ii,1) = I1;
  end
  if ~isempty(I2)
    Igold(ii,2) = I2;
  end
end

figure
for ii = 1:5
  subplot(2,3,ii),hold on
  plot(Chromatograms(Igold(ii,1),:))
  plot(Chromatograms(Igold(ii,2),:))
  title(num2str(scoreMatrix(Igold(ii,1),Igold(ii,2))))
  legend(nameGold{ii,1},nameGold{ii,2})
  xlim([0 66])
end
%subplot(2,3,1),title('Gold standards')

maxscores = sort(scoreMatrix(:),'descend');
maxscores = maxscores([1 5 10 15 20 25]);
figure
for ii = 1:6
  I = find(scoreMatrix(:) == maxscores(ii));
  [I1,I2] = ind2sub(size(scoreMatrix),I);
  I1 = I1(randsample(length(I1),1));
  I2 = I2(randsample(length(I2),1));
  subplot(2,3,ii),hold on
  plot(Chromatograms(I1,:)')
  plot(Chromatograms(I2,:)')
  title(num2str(scoreMatrix(I1,I2)))
  legend(tmp{I1},tmp{I2})
end

figure
I1 = possibleInts(:) & TP_Matrix(:)==1;
I2 = possibleInts(:) & FP_Matrix(:)==1;
fn = fieldnames(Dist);
for ii = 1:length(fn)
  subplot(3,3,ii),hold on
  x1 = min([Dist.(fn{ii})(I1)' Dist.(fn{ii})(I2)']);
  x2 = max([Dist.(fn{ii})(I1)' Dist.(fn{ii})(I2)']);
  xi = linspace(x1,x2,101);
  if ii == 4
    xi = 1:8;
  end
  h1 = hist(Dist.(fn{ii})(I1),xi);
  h2 = hist(Dist.(fn{ii})(I2),xi);
  plot(xi,h1/sum(h1))
  plot(xi,h2/sum(h2))
  xlim([x1 x2])
  title(fn{ii})
  
  if ii==1
    legend('Pos.','Neg.')
  end
end



%% Make PR curve from different Dist fields

figure

fn = fieldnames(Dist);
I1 = possibleInts(:) & TP_Matrix(:)==1;
I2 = possibleInts(:) & FP_Matrix(:)==1;
for ii = 1:length(fn)
  
  x1 = min([Dist.(fn{ii})(I1)' Dist.(fn{ii})(I2)']);
  x2 = max([Dist.(fn{ii})(I1)' Dist.(fn{ii})(I2)']);
  xi = linspace(x1,x2,101);
  
  prec = nan(size(xi));
  rec = nan(size(xi));
  Nint = nan(size(xi));
  for jj = 1:length(xi)
    TP = sum(Dist.(fn{ii})(:) < xi(jj) & I1);
    FP = sum(Dist.(fn{ii})(:) < xi(jj) & I2);
    FN = sum(Dist.(fn{ii})(:) >= xi(jj) & I1);
    prec(jj) = TP / (TP + FP);
    rec(jj) = TP / (TP + FN);
    Nint(jj) = TP + FP;
  end
  
  subplot(3,3,ii)
  plot(rec,prec)
  axis([0 1 0 1])
  title(fn{ii})
end



%%  Mike Carlson's myterious data: Figure out why Dist.(all fields) does worse than just Dist.R2
% Look at feats.

ff = '/Users/Mercy/Academics/Foster/MikeCarlson/Data_files/feats_20160710.mat';
load(ff)
fn{11} = 'All';

figure
bar([Ninter Ninter_all])
set(gca,'xtick',1:length(Ninter)+1,'xticklabel',fn)

feats_individually = zeros(15,10);
for ii = 1:10
  feats_individually(:,ii) = Feats{ii};
end

s1 = std(Feats_all) / sqrt(15);
m1 = mean(Feats_all);
I = repmat((1:10),15,1);

figure,hold on
%errorbar(1:10,m2,m2-s2,m2+s2)
errorbar(1:10,m1,m1-s1,m1+s1,'m')
boxplot(feats_individually(:),I(:))
legend('Together')
set(gca,'xtick',1:length(Ninter)+1,'xticklabel',fn)


figure,subplot(1,2,1)
imagesc(Feats_all>5)
colorbar;
cx = caxis;
set(gca,'xtick',1:length(Ninter)+1,'xticklabel',fn);
title('Together'),
subplot(1,2,2),
imagesc(feats_individually>5),
caxis(cx);
colorbar,
title('By themselves'),
set(gca,'xtick',1:length(Ninter)+1,'xticklabel',fn)



%%  Mike Carlson's myterious data: play around with pairwise file
% ROC_PCPSILAC has only Dist.R2 selected

stringRange = [300 500 700 800 825 850 900 950 975 990 995 997 998];
Ninter = [129 136 149 251 221 252 251 210 149 150 130 130 136];

figure
plot(stringRange,Ninter)
xlabel('STRING cutoff')
ylabel('Number of predicted interactions @ precision = 50%')


% tissue data
Ninteractions_all = [13977 21681 20273 16741 16827 24449];
Ninteractions_justR2 = [14057 20229 20293 24080 18549 25158];



%%  Mike Carlson's myterious data: Possible explanation
% Two things going on:
%   1. R^2 does better by itself. This may be a feature selection problem (asked on SE). Ignore now
%   2. STRING is not as informative as CORUM
%
% Compare string-separation to corum-separation
% For each feature, calculate AUC between feature(class==1) and feature(class==0)

% CORUM-separation for tissue data
dir1 = '/Users/Mercy/Academics/Foster/Tissue_PCPSILAC/PCPSILAC_analysis/';
featauc = zeros(6,9);
dm = zeros(size(featauc));
for ii = 1:6
  ii
  sf = [dir1 '/Output/tmp/' 'score_rep' num2str(ii) '.mat'];
  load(sf)
  I = possibleInts(:);
  fn = fieldnames(Dist);
  for jj = 1:length(fn)
    jj
    clear x
    x(:,1) = Dist.(fn{jj})(I);
    x(:,2) = TP_Matrix(I)==1;
    I2 = ~isnan(x(:,1)) & ~isinf(x(:,1));
    x = x(I2,:);
    tmp = roc(x);
    featauc(ii,jj) = tmp.AUC;
    
    dm(ii,jj) = nanmean(x(x(:,2)==0,1)) - nanmean(x(x(:,2)==1,1));
  end
end


dir2 = '/Users/Mercy/Academics/Foster/MikeCarlson/PCPSILAC_analysis/';
sf = [dir2 '/Output/tmp/' 'score_rep1.mat'];
load(sf)
I = possibleInts(:);
fn = fieldnames(Dist);
featauc_string = zeros(9,1);
dm_string = zeros(9,1);
for jj = 1:length(fn)
  jj
  clear x
  x(:,1) = Dist.(fn{jj})(I);
  x(:,2) = TP_Matrix(I)==1;
  I2 = ~isnan(x(:,1)) & ~isinf(x(:,1));
  x = x(I2,:);
  x = x(1:10:end,:);
  tmp = roc(x);
  featauc_string(jj) = tmp.AUC;
  
  dm_string(jj) = nanmean(x(x(:,2)==0,1)) - nanmean(x(x(:,2)==1,1));
end


figure,

subplot(2,1,1),hold on
[~,I] = sort(featauc_string,'descend');
X = repmat(1:9,6,1);
X = X(:,:);
y = featauc(:,I);
scatter(X(:),y(:))
plot(featauc_string(I))
set(gca,'xtick',1:9,'xticklabel',fn(I))
legend('CORUM','STRING')

subplot(2,1,2),hold on
[~,I] = sort(dm_string,'descend');
X = repmat(1:9,6,1);
X = X(:,:);
y = dm(:,I);
scatter(X(:),y(:))
plot(dm_string(I))
set(gca,'xtick',1:9,'xticklabel',fn(I))
legend('CORUM','STRING')



%% Woo! Try a decoy search to estimate FPs

sf_real = '/Users/Mercy/Academics/Foster/Tissue_PCPSILAC/PCPSILAC_analysis//Output/tmp/score_rep1.mat';
load(sf_real)
score_real = scoreMatrix;

sf_shuffled = '/Users/Mercy/Academics/Foster/Tissue_PCPSILAC/PCPSILAC_analysis//Output/tmp/score_rep1_shuffledIDs.mat';
load(sf_shuffled)
score_shuffled = scoreMatrix;

ds = linspace(0,1,501);
prec = zeros(size(ds));
Ninter = zeros(size(ds));
for ii = 1:length(ds)
  ii
  FP = sum(score_shuffled > ds(ii));
  Ninter(ii) = sum(score_shuffled > ds(ii) | score_real > ds(ii));
  prec(ii) = 1 - FP / Ninter(ii);
end

x = linspace(0,1,501);
h1 = hist(score_real,x);
h2 = hist(score_shuffled,x);
figure,hold on
plot(x,h1)
plot(x,h2)
legend('Real','Shuffled')




%% Reference-less FPs, take 2

% Make R_raw and p_raw for any replicate

I = Rraw==1 | Rraw==-1;
Rraw(I) = nan;
praw(I) = nan;
N = sum(~isnan(Chromatograms_raw),2);
I = N<=10;
Rraw(I,I) = nan;
praw(I,I) = nan;
Rraw(abs(Rraw)>.999) = nan;

ds = exp(-(0:50));
Ninter = zeros(size(ds));
fdr = zeros(size(ds));
for ii = 1:length(ds)
  I = praw(:)<ds(ii);
  TP = sum(I & Rraw(:)>0);
  FP = sum(I & Rraw(:)<0);
  
  Ninter(ii) = TP + FP;
  fdr(ii) = FP / (TP + FP);
end


figure
plot(fdr,Ninter)


[~,I] = min(abs(fdr - 0.01));
pcut01 = ds(I);
[~,I] = min(abs(fdr - 0.05));
pcut05 = ds(I);

I01 = praw(:)<pcut01;
I05 = praw(:)>=pcut01 & praw(:)<pcut05;

x = linspace(-1.1,1.1,101);
h01 = hist(Rraw(I01),x);
h05 = hist(Rraw(I05),x);
h0 = hist(Rraw(~I01 & ~I05),x);

figure
b = bar(x,[h01;h05;h0]','stacked');
b(1).FaceColor = 'c';
b(1).EdgeColor = 'c';
b(2).FaceColor = 'b';
b(2).EdgeColor = 'b';
b(3).FaceColor = [.7 .7 .7];
b(3).EdgeColor = [.7 .7 .7];
xlabel('Correlation coefficient','fontsize',14)
ylabel('Count','fontsize',14)
xlim([-1 1])
s1 = ['1% FDR, N = ' num2str(sum(I01)/2)];
s2 = ['5% FDR, N = ' num2str(sum(I05 | I01)/2)];
legend(s1,s2,'Non-sig.')
set(gca,'fontsize',14)



%% Now, compare the predicted interactors to CORUM

jj = 3;
clear Chromatograms_raw Dist TP_Matrix FP_Matrix possibleInts
%sf = ['/Users/Mercy/Academics/Foster/Tissue_PCPSILAC/PCPSILAC_analysis//Output/tmp/score_rep' num2str(jj) '.mat'];
%sf = ['/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/Output/tmp/score_rep' num2str(jj) '.mat'];
sf = '/Users/Mercy/Academics/Foster/MikeCarlson/PCPSILAC_analysis//Output/tmp/score_rep1.mat';
load(sf)
%RR = reshape(RR,sqrt(length(RR)),sqrt(length(RR)));

%[R,p] = corrcoef(Chromatograms_raw','rows','pairwise');

% Remove nans from Rpraw, since can never be found to interact
% I = ~isnan(Dist.Rpraw);
% Dist.Rpraw = Dist.Rpraw(I);
% TP_Matrix = TP_Matrix(I);
% possibleInts = possibleInts(I);
% FP_Matrix = FP_Matrix(I);
% RR = RR(I);

randTP_Matrix = zeros(size(TP_Matrix));
I = randsample(length(TP_Matrix(:)==1 & possibleInts(:)==1),sum(TP_Matrix(:)==1 & possibleInts(:)==1)*1);
randTP_Matrix(I) = 1;

ds = exp(-(0:2:200));
recallTP_cor = zeros(size(ds));
recallFP_cor = zeros(size(ds));
recallFPrand = zeros(size(ds));
fdr_apms = zeros(size(ds));
prec_cor = zeros(size(ds));
Ninter = zeros(size(ds));
for ii = 1:length(ds)
  ii
  % CORUM TP
  I_data = Dist.Rpraw<ds(ii);
  I_corum = TP_Matrix==1 & possibleInts==1;
  I_overlap = I_data & I_corum;
  recallTP_cor(ii) = sum(I_overlap(:)) / sum(I_corum(:));
  
  % fdr, CORUM method
  TP = sum(Dist.Rpraw(:)<ds(ii) & TP_Matrix(:)==1 & possibleInts(:)==1);
  FP = sum(Dist.Rpraw(:)<ds(ii) & TP_Matrix(:)==0 & possibleInts(:)==1);
  prec_cor(ii) = TP / (TP + FP);
  
  % CORUM "FP"
  I_corum = FP_Matrix==1 & possibleInts==1;
  I_overlap = I_data & I_corum;
  recallFP_cor(ii) = sum(I_overlap(:)) / sum(I_corum(:));
  
  % Random "FP"
  I_overlap = I_data & randTP_Matrix==1;
  recallFPrand(ii) = sum(I_overlap(:)) / sum(randTP_Matrix(:));
  
  % fdr, AP-MS method
  TP = sum(Dist.Rpraw(:)<ds(ii) & RR(:)>0);
  FP = sum(Dist.Rpraw(:)<ds(ii) & RR(:)<0);
  fdr_apms(ii) = FP / (TP + FP);
  
  Ninter(ii) = TP + FP;
end

figure
subplot(2,2,1)
semilogx(ds,recallTP_cor)
hold on
semilogx(ds,recallFP_cor,'r')
semilogx(ds,recallFPrand,'m')
%plot([1 1]*(pcut01),y,'--k')
%plot([1 1]*(pcut05),y,'--k')
legend('CORUM TP','CORUM FP','Random','location','northwest')
grid on
set(gca,'fontsize',14)
xlabel('P-value threshold','fontsize',14)
ylabel('Recall','fontsize',14)
x = xlim;

subplot(2,2,3)
semilogx(ds,1-prec_cor)
%hold on
%semilogx(ds,1-prec_cor)
xlim(x)
grid on
legend('CORUM TP','location','northwest')
ylabel('FDR','fontsize',14)
xlabel('P-value threshold','fontsize',14)


% Make precision-recall curve from score
I = possibleInts(:)==1;
score = scoreMatrix(I);
class = TP_Matrix(I);
[~, tmp] = calcPPIthreshold(score,class,0.9);

subplot(2,2,[2 4]),hold on
plot(recallTP_cor,prec_cor)
plot(tmp.recRange,tmp.precRange,'r')
title('Precision-Recall, CORUM TP')
legend('AP-MS method','Our method')
ylabel('Precision','fontsize',14)
xlabel('Recall','fontsize',14)
mySound
pause(.01)


%% Try to explain the hump at R = -0.5

x = 1:55;
rr = nan(5000,1);
for iter = 1:5000
  y1 = normpdf(x,ceil(rand*55),ceil(rand*10));
  y1 = y1 / max(y1);
  I = randsample(55,5);
  y1(I) = nan;
  y2 = normpdf(x,ceil(rand*55),ceil(rand*10));
  y2 = y2 / max(y2);
  I = randsample(55,5);
  y2(I) = nan;
  y1(y1<.1) = nan;
  y2(y2<.1) = nan;
  y1 = y1+rand(size(y1))*.95;
  y2 = y2+rand(size(y2))*.95;
  I = ~isnan(y1) & ~isnan(y2);
  if sum(I) ==0
    continue
  end
  rr(iter) = corr(y1(I)',y2(I)');
end

figure,hist(rr,linspace(-1,1,26))



%% Predict Mike Carlson's interactions with AP-MS method

sf = '/Users/Mercy/Academics/Foster/MikeCarlson/PCPSILAC_analysis//Output/tmp/score_rep1.mat';



%% Sanity check: compare all-replicate-scoreMatrix to just R^2

% scoreMatrix
%[xcutoff, tmp] = calcPPIthreshold(scoreMatrix(possList==1),classList(possList==1)',desiredPrecision);

N = 51;
dR = 10.^(linspace(-5,0,N));
rec1 = nan(N,1);
prec1 = nan(N,1);
rec2 = nan(N,1);
prec2 = nan(N,1);
rec3 = nan(N,1);
prec3 = nan(N,1);
dS = 1 - 40.^(linspace(-5,0,N));
rec0 = nan(N,1);
prec0 = nan(N,1);
for ii = 1:length(dR)
  ii
  TP = sum(possList'==1 & classList'==1 & DistList(:,2)<dR(ii));
  FP = sum(possList'==1 & classList'==0 & DistList(:,2)<dR(ii));
  FN = sum(possList'==1 & classList'==1 & DistList(:,2)>dR(ii));
  prec1(ii) = TP / (TP + FP);
  rec1(ii) = TP / (TP + FN);
  
  TP = sum(possList'==1 & classList'==1 & DistList(:,7)<dR(ii));
  FP = sum(possList'==1 & classList'==0 & DistList(:,7)<dR(ii));
  FN = sum(possList'==1 & classList'==1 & DistList(:,7)>dR(ii));
  prec2(ii) = TP / (TP + FP);
  rec2(ii) = TP / (TP + FN);
  
  TP = sum(possList'==1 & classList'==1 & DistList(:,12)<dR(ii));
  FP = sum(possList'==1 & classList'==0 & DistList(:,12)<dR(ii));
  FN = sum(possList'==1 & classList'==1 & DistList(:,12)>dR(ii));
  prec3(ii) = TP / (TP + FP);
  rec3(ii) = TP / (TP + FN);
  
  TP = sum(possList'==1 & classList'==1 & score>dS(ii));
  FP = sum(possList'==1 & classList'==0 & score>dS(ii));
  FN = sum(possList'==1 & classList'==1 & score<dS(ii));
  prec0(ii) = TP / (TP + FP);
  rec0(ii) = TP / (TP + FN);
end

figure,hold on,
plot(rec0,prec0,'linewidth',2,'color','k')
plot(rec1,prec1,'g')
plot(rec2,prec2,'r')
plot(rec3,prec3,'m')
legend('classifier','R2 1','R2 2 ','R2 3')


%% Something is off
% test the idea that NEW interactions for each replicate are good
% but the overlapping ones are bad

Ifirstnan = [1 0 0 0];
for ii = 1:3
  Ifirstnan(ii+1) = find(~isnan(DistList(:,ii*5 -3)),1,'last');
end

I = 1:size(DistList,1);
dR = 10.^(linspace(-5,0,51));

for ii = 1:3
  I0 = Ifirstnan(ii) : Ifirstnan(ii+1);
  I1 = ismember(I,I0); % new
  
  figure,hold on
  rec1 = zeros(size(dR));
  rec2 = zeros(size(dR));
  prec1 = zeros(size(dR));
  prec2 = zeros(size(dR));
  for jj = 1:length(dR)
    jj
    
    % New interactions
    TP = nansum(possList(I1)'==1 & classList(I1)'==1 & DistList(I1,2 + (ii-1)*5)<dR(jj));
    FP = nansum(possList(I1)'==1 & classList(I1)'==0 & DistList(I1,2 + (ii-1)*5)<dR(jj));
    FN = nansum(possList(I1)'==1 & classList(I1)'==1 & DistList(I1,2 + (ii-1)*5)>dR(jj));
    prec1(jj) = TP / (TP + FP);
    rec1(jj) = TP / (TP + FN);
    
    % Overlapping interactions
    TP = nansum(possList(~I1)'==1 & classList(~I1)'==1 & DistList(~I1,2 + (ii-1)*5)<dR(jj));
    FP = nansum(possList(~I1)'==1 & classList(~I1)'==0 & DistList(~I1,2 + (ii-1)*5)<dR(jj));
    FN = nansum(possList(~I1)'==1 & classList(~I1)'==1 & DistList(~I1,2 + (ii-1)*5)>dR(jj));
    prec2(jj) = TP / (TP + FP);
    rec2(jj) = TP / (TP + FN);
  end
  
  plot(rec1,prec1,'g')
  plot(rec2,prec2,'r')
  legend('new','overlapping')
  title(num2str(replicatesThisChannel(ii)))
  grid on
  pause(.001)
end


%% Why is Inew worse than ~Inew??

Ifirstnan = [1 0 0 0];
for ii = 1:3
  Ifirstnan(ii+1) = find(~isnan(DistList(:,ii*5 -3)),1,'last') + 1;
end

ii = 2;
I = 1:size(DistList,1);
I0 = Ifirstnan(ii) : Ifirstnan(ii+1);
I1 = ismember(I,I0); % new
  
figure
subplot(2,1,1)
hist(DistList(classList==1 & ~I1,7))
title('Positive class, overlapping')
subplot(2,1,2)
hist(DistList(classList==1 & I1,7))
title('Positive class, new')


%% Manual check
% get Inew for replicate 6, replicatesThisChannel = [5 6 4]
% keep [indexList(I1,:) classList(I1)' possList(I1)' DistList(I1,:)]
% compare to the same interactions, but for [6 5 4], replicate 6


% * run w/ replicatesThisChannel = [4 5 6]
% Matnew5 = [indexList(I1,:) classList(I1)' possList(I1)' DistList(I1,:)];
% DistList5 = DistList;
% protAll_save = proteinAll;

% * run w/ replicatesThisChannel = [6 5 4]

% try it for a single line
I = find(Matnew5(:,4)==1)';
Imat = zeros(size(I));
cc = 0;
for ii = I
  ii
  cc = cc+1;
  i1 = Matnew5(ii,1);
  i2 = Matnew5(ii,2);
  protA = protAll_save{i1};
  protB = protAll_save{i2};
  i3(1) = find(ismember(proteinAll,protA));
  i3(2) = find(ismember(proteinAll,protB));
  i3 = sort(i3);
  Imat(cc) = find(indexList(:,1)==i3(1) & indexList(:,2)==i3(2));
  %tmp1 = Matnew5(ii,[3 4 10:14]);
  %tmp2 = [classList(Imat) possList(Imat) DistList(Imat,6:10)];
  
  %if sum(tmp1 - tmp2)~=0
  %  pause
  %end
end



%% Oooookay.
% I'm making DistList correctly.
% But then why does the R^2 outperform the score?
% It must be a badly performing classifier

% try it without any nans
%I = find(sum(isnan(DistList),2) == 0);
I = 1:size(DistList,1);

score2 = scorenb(DistList(I,:),possList(I),classList(I));
score2 = nanmedian(score2,2);


N = 51;
dR = 10.^(linspace(-5,0,N));
dS = 1 - 100.^(linspace(-5,0,N));
clear rec1 prec1 rec2 prec2 rec3 prec3
for ii = 1:length(dR)
  ii
  TP = sum(possList(I)'==1 & classList(I)'==1 & DistList(I,2)<dR(ii));
  FP = sum(possList(I)'==1 & classList(I)'==0 & DistList(I,2)<dR(ii));
  FN = sum(possList(I)'==1 & classList(I)'==1 & DistList(I,2)>dR(ii));
  prec1(ii) = TP / (TP + FP);
  rec1(ii) = TP / (TP + FN);
  
  TP = sum(possList(I)'==1 & classList(I)'==1 & DistList(I,7)<dR(ii));
  FP = sum(possList(I)'==1 & classList(I)'==0 & DistList(I,7)<dR(ii));
  FN = sum(possList(I)'==1 & classList(I)'==1 & DistList(I,7)>dR(ii));
  prec2(ii) = TP / (TP + FP);
  rec2(ii) = TP / (TP + FN);
  
  TP = sum(possList(I)'==1 & classList(I)'==1 & DistList(I,12)<dR(ii));
  FP = sum(possList(I)'==1 & classList(I)'==0 & DistList(I,12)<dR(ii));
  FN = sum(possList(I)'==1 & classList(I)'==1 & DistList(I,12)>dR(ii));
  prec3(ii) = TP / (TP + FP);
  rec3(ii) = TP / (TP + FN);
  
  TP = sum(possList(I)'==1 & classList(I)'==1 & score2>dS(ii));
  FP = sum(possList(I)'==1 & classList(I)'==0 & score2>dS(ii));
  FN = sum(possList(I)'==1 & classList(I)'==1 & score2<dS(ii));
  prec0(ii) = TP / (TP + FP);
  rec0(ii) = TP / (TP + FN);
end

figure,hold on,
plot(rec0,prec0,'linewidth',2,'color','k')
plot(rec1,prec1,'g')
plot(rec2,prec2,'r')
plot(rec3,prec3,'m')
legend('classifier','R2 1','R2 2 ','R2 3')
grid on


%% Check Ali's skewed histograms

% make ratio data with mean 2
yy = ones(10000,1);
xx = yy + randn(10000,1);

data = exp(xx./yy)*2;
log2data = log2(data);
dataNorm = log2data - nanmedian(log2data);

figure
subplot(1,3,1)
hist(data,0:.5:300)
xlim([-1 30])
title('data')
subplot(1,3,2)
hist(log2data,-5:.5:8)
xlim([-5 8])
title('log2 data')
subplot(1,3,3)
hist(dataNorm,-5:.5:8)
xlim([-5 8])
title('log2 normalized data')


%% Check the output of Anna's file

fn = '/Users/Mercy/Academics/Foster/Anna/shRNA searching/AW_TRC_9K_locations.csv';
fid = fopen(fn,'r');
fgetl(fid);
cc = 0;
ids9k = cell(9656,1);
nn9k = zeros(9656,1);
while ~feof(fid)
  t = fgetl(fid);
  t1 = strsplit(t,',');
  cc = cc+1;
  ids9k{cc} = t1{1};
  nn9k(cc) = length(t1);
end
ids9k = ids(1:cc);
nn9k = nn9k(1:cc);


fn = '/Users/Mercy/Academics/Foster/Anna/shRNA searching/AW_id_locations.csv';
fid = fopen(fn,'r');
fgetl(fid);
cc = 0;
ids3k = cell(9656,1);
nn3k = zeros(9656,1);
while ~feof(fid)
  t = fgetl(fid);
  t1 = strsplit(t,',');
  cc = cc+1;
  ids3k{cc} = t1{1};
  nn3k(cc) = length(t1);
end
ids3k = ids3k(1:cc);
nn3k = nn3k(1:cc);