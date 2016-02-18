% Figures to make:
% Fig. 1a: Fold change. Venn diagram of proteins where a change was detected by i) Nick and ii) me.
%   How much overlap is there?
% Fig. 1b: Fold change. Scatter plot of fold-change-Nick vs fold-change-Greg.
% Fig. 2a: Ttest / ranksum. Similar to 1a except for the method of change detection.


%% Fold change

% 1. Read in proteins where Nick detected a change.
fn = '/Users/Mercy/Academics/Foster/NickCodeData/3_Comparsion processing/Analysis output/Unique_gaussians_with_changes.csv';
data = importdata(fn);
Nprot = size(data.data,1) - 1;
Protein_name = cell(Nprot,1);
Change = zeros(Nprot,1);
Center = zeros(Nprot,1);
Height = zeros(Nprot,1);
Fold = zeros(Nprot,1);
for ri = 2:Nprot+1
  tmp = data.textdata{ri,8};
  if strcmp(tmp,'No change')
    Change(ri-1) = 0;
  elseif strcmp(tmp,'Decrease')
    Change(ri-1) = -1;
  elseif strcmp(tmp,'Increase')
    Change(ri-1) = 1;
  else
    Change(ri-1) = nan;
  end
  Protein_name{ri-1} = data.textdata{ri,1};
  Center(ri) = str2num(data.textdata{ri,3});
  Height(ri) = str2num(data.textdata{ri,4});
  Fold(ri) = data.data(ri,3);
end

% 2. Read in proteins where I detected a change
fn = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/Data/Comparison/Unique_gaussians_with_changes_betweenMvsLandHvslb.csv';
data = importdata(fn);
Nprot = size(data.data,1) - 1;
Protein_nameg = cell(Nprot,1);
Changeg = zeros(Nprot,1);
Centerg = zeros(Nprot,1);
Heightg = zeros(Nprot,1);
Foldg = zeros(Nprot,1);
for ri = 2:Nprot+1
  tmp = data.textdata{ri,8};
  if strcmp(tmp,'No change')
    Changeg(ri-1) = 0;
  elseif strcmp(tmp,'Decrease')
    Changeg(ri-1) = -1;
  elseif strcmp(tmp,'Increase')
    Changeg(ri-1) = 1;
  else
    Changeg(ri-1) = nan;
  end
  Protein_nameg{ri-1} = data.textdata{ri,1};
  Centerg(ri) = str2num(data.textdata{ri,3});
  Heightg(ri) = str2num(data.textdata{ri,2});
  Foldg(ri) = data.data(ri,1);
end

% 3. Make that Venn diagram.
% for each of Nick's Gaussians
%   find the closest of my Gaussians w/ DeltaCenter<3
foldcheck = zeros(length(Protein_name),8);
for ri = 1:length(Protein_name)
  I = find(ismember(Protein_nameg,Protein_name{ri}));
  dC = abs(Centerg(I) - Center(ri));
  [I2,I3] = min(dC);
  if I2<2
    I4 = I(I3);
    foldcheck(ri,:) = [Change(ri) Changeg(I4) Fold(ri) Foldg(I4) Center(ri) Height(ri) Centerg(I4) Heightg(I4)];
  else
    foldcheck(ri,:) = nan(8,1);
  end
end

Incn = find(foldcheck(:,1)==1);
Incg = find(foldcheck(:,2)==1);
Inc2 = find(foldcheck(:,1)==1 & foldcheck(:,2)==1);
Decn = find(foldcheck(:,1)==-1);
Decg = find(foldcheck(:,2)==-1);
Dec2 = find(foldcheck(:,1)==-1 & foldcheck(:,2)==-1);
NCn = find(foldcheck(:,1)==0);
NCg = find(foldcheck(:,2)==0);
NC2 = find(foldcheck(:,1)==0 & foldcheck(:,2)==0);


figure
subplot(3,1,1)
myVenn2([length(Incn) length(Incg)], length(Inc2))
text(0,0,num2str(length(Inc2)))
text(-15,10,num2str(length(Incn) - length(Inc2)))
text(16,10,num2str(length(Incg) - length(Inc2)))
set(gca,'xtick',[],'ytick',[])
title({'Change detected with fold>2' 'Increase'})
subplot(3,1,2)
myVenn2([length(NCn) length(NCg)], length(NC2))
text(0,0,num2str(length(NC2)))
text(-50,37,num2str(length(NCn) - length(NC2)))
text(50,35,num2str(length(NCg) - length(NC2)))
set(gca,'xtick',[],'ytick',[])
title('No change')
subplot(3,1,3)
myVenn2([length(Decn) length(Decg)], length(Dec2))
text(0,0,num2str(length(Dec2)))
text(-10,10,num2str(length(Decn) - length(Dec2)))
text(10,10,num2str(length(Decg) - length(Dec2)))
set(gca,'xtick',[],'ytick',[])
title('Decrease')
set(gcf,'units','normalized','position',[.1 .1 .3 .8])


figure
subplot(2,1,1),hold on
I = (foldcheck(:,3)>1 & foldcheck(:,4)>1) | (foldcheck(:,3)<-1 & foldcheck(:,4)<-1) | abs(foldcheck(:,3))<1 & abs(foldcheck(:,4))<1;
I2 = ~I & ~isnan(foldcheck(:,1));
scatter(foldcheck(I,3),foldcheck(I,4),6,'b','filled')
scatter(foldcheck(I2,3),foldcheck(I2,4),6,'r','filled')
text(-2.7,4.5,['N = ' num2str(sum(~isnan(foldcheck(:,1)))) ' Gaussians'])
legend(['Same, N=' num2str(sum(I))],['Diff, N=' num2str(sum(I2))],'location','southeast')
xlabel('Fold change, old code','fontsize',14)
ylabel('Fold change, new code','fontsize',14)
title('Fold Change, Unique_Gaussians_with_changes.csv')
grid on
subplot(2,1,2),hold on
values = log10(hist3([foldcheck(:,3) foldcheck(:,4)],{linspace(-3,5,201),linspace(-3,5,201)}));
imagesc(linspace(-3,5,201),linspace(-3,5,201),values)
axis xy
axis([-3 5 -3 5])
plot([-3 5],[-3 5],'--r')
xlabel('Fold change, old code','fontsize',14)
ylabel('Fold change, new code','fontsize',14)
title('Log Density')
set(gcf,'units','normalized','position',[.1 .1 .3 .8])



%% Two example proteins

badProts = (~I' & ~isnan(foldcheck(:,1))' & abs(diff(foldcheck(:,3:4)'))>1);



%% T-test, ranksum

% 1. Load Nick's data
fn = '/Users/Mercy/Academics/Foster/NickCodeData/3_Comparsion processing/Analysis output/Perseus_enrichment_Gaussian_level_file.csv';
fid = fopen(fn);
line = fgetl(fid); %header
Center = zeros(5000,1);
Protein_name = cell(5000,1);
ttp = nan(5000,1);
wwp = nan(5000,1);
ttChange = zeros(5000,1);
wwChange = zeros(5000,1);
foldChange = zeros(5000,1);
jj = 0;
while ~feof(fid)
  jj = jj + 1;
  line = strsplit(fgetl(fid),',');
  Center(jj) = str2num(line{1});
  Protein_name{jj} = line{2};
  foldChange(jj) = str2num(line{3});
  try ttp(jj) = str2num(line{11});end
  try wwp(jj) = str2num(line{13});end
  wwChange(jj) = 0;
  ttChange(jj) = 0;
  if length(line)>11
    if strcmp('+',line{12})
      ttChange(jj) = 1;
    end
  end
  if length(line)>13
    if strcmp('+',line{14})
      wwChange(jj) = 1;
    end
  end
end
fclose(fid);


% 2. Load My data
fn = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/Data/Comparison/Perseus_enrichment_Gaussian_level_fileb.csv';
fid = fopen(fn);
line = fgetl(fid); %header
Centerg = zeros(5000,1);
Protein_nameg = cell(5000,1);
ttpg = nan(5000,1);
wwpg = nan(5000,1);
ttChangeg = zeros(5000,1);
wwChangeg = zeros(5000,1);
foldChangeg = zeros(5000,1);
jj = 0;
while ~feof(fid)
  jj = jj + 1;
  line = strsplit(fgetl(fid),',');
  Centerg(jj) = str2num(line{1});
  Protein_nameg{jj} = line{2};
  foldChangeg(jj) = str2num(line{3});
  try ttpg(jj) = str2num(line{6});end
  try wwpg(jj) = str2num(line{8});end
  wwChangeg(jj) = 0;
  ttChangeg(jj) = 0;
  if length(line)>11
    if strcmp('+',line{12})
      ttChangeg(jj) = 1;
    end
  end
  if length(line)>13
    if strcmp('+',line{14})
      wwChangeg(jj) = 1;
    end
  end
end
fclose(fid);


% 3. Make that Venn diagram.
% for each of Nick's Gaussians
%   find the closest of my Gaussians w/ DeltaCenter<3
pcheck = zeros(length(Protein_name),8);
foldcheck2 = zeros(length(Protein_name),2);
for ri = 1:length(Protein_name)
  I = find(ismember(Protein_nameg,Protein_name{ri}));
  dC = abs(Centerg(I) - Center(ri));
  [I2,I3] = min(dC);
  if I2<2
    I4 = I(I3);
    pcheck(ri,:) = [ttChange(ri) ttChangeg(I4) wwChange(ri) wwChangeg(I4) ttp(ri) ttpg(I4) wwp(ri) wwpg(I4) ];
    foldcheck2(ri,:) = [foldChangeg(I4) foldChange(ri)];
  else
    pcheck(ri,:) = nan(8,1);
  end
end

ttn = find(pcheck(:,1)==1);
ttg = find(pcheck(:,2)==1);
tt2 = find(pcheck(:,1)==1 & pcheck(:,2)==1);

wwn = find(pcheck(:,3)==1);
wwg = find(pcheck(:,4)==1);
ww2 = find(pcheck(:,3)==1 & pcheck(:,4)==1);

ttpn = find(pcheck(:,5)==1);
ttpg = find(pcheck(:,6)==1);
ttp2 = find(pcheck(:,5)==1 & pcheck(:,6)==1);

wwpn = find(pcheck(:,7)==1);
wwpg = find(pcheck(:,8)==1);
wwp2 = find(pcheck(:,7)==1 & pcheck(:,8)==1);


figure
subplot(2,1,1)
myVenn2([length(ttn) length(ttg)], length(tt2))
text(0,0,num2str(length(ttn)))
text(-8,10,num2str(length(ttn) - length(tt2)))
text(11,9,num2str(length(ttg) - length(tt2)))
set(gca,'xtick',[],'ytick',[])
title('Change detected with ttest')
subplot(2,1,2)
myVenn2([length(wwn) length(wwg)], length(ww2))
text(0,0,num2str(length(wwn)))
text(-5,5.2,num2str(length(wwn) - length(ww2)))
text(7,5,num2str(length(wwg) - length(ww2)))
set(gca,'xtick',[],'ytick',[])
title('Change detected with rank-sum')
set(gcf,'units','normalized','position',[.1 .1 .3 .8])

figure
subplot(2,2,1),hold on
I = ((pcheck(:,1) & pcheck(:,2)) | (~pcheck(:,1) & ~pcheck(:,2)));
I2 = ~I & ~isnan(pcheck(:,1));
scatter(log(pcheck(I,5)),log(pcheck(I,6)),6,'b','filled')
scatter(log(pcheck(I2,5)),log(pcheck(I2,6)),6,'r','filled')
ax1=axis;
xlabel('T-test log p-value, old code','fontsize',14)
ylabel('T-test log p-value, new code','fontsize',14)
legend(['Same, N=' num2str(sum(I))],['Diff, N=' num2str(sum(I2))],'location','southeast')
title('Scatter')
grid on
axis(ax1)
subplot(2,2,2),hold on
I = ((pcheck(:,3) & pcheck(:,4)) | (~pcheck(:,3) & ~pcheck(:,4)));
I2 = ~I & ~isnan(pcheck(:,1));
scatter(log(pcheck(I,7)),log(pcheck(I,8)),6,'b','filled')
scatter(log(pcheck(I2,7)),log(pcheck(I2,8)),6,'r','filled')
ax2=axis;
xlabel('Mann-Whitney log p-value, old code','fontsize',14)
ylabel('Mann-Whitney log p-value, new code','fontsize',14)
legend(['Same, N=' num2str(sum(I))],['Diff, N=' num2str(sum(I2))],'location','southeast')
title('Scatter')
grid on
axis(ax2)
subplot(2,2,3),hold on
values = log10(hist3(log([pcheck(:,5) pcheck(:,6)]),{linspace(ax1(1),ax1(2),101),linspace(ax1(1),ax1(2),101)}));
imagesc(linspace(ax1(1),ax1(2),101),linspace(ax1(1),ax1(2),101),values)
axis xy
axis(ax1)
xlabel('T-test log p-value, old code','fontsize',14)
ylabel('T-test log p-value, new code','fontsize',14)
title('Log Density')
subplot(2,2,4),hold on
values = log10(hist3(log([pcheck(:,7) pcheck(:,8)]),{linspace(ax2(1),ax2(2),101),linspace(ax2(1),ax2(2),101)}));
imagesc(linspace(ax2(1),ax2(2),101),linspace(ax2(1),ax2(2),101),values)
axis xy
axis(ax2)
xlabel('T-test log p-value, old code','fontsize',14)
ylabel('T-test log p-value, new code','fontsize',14)
title('Log Density')
set(gcf,'units','normalized','position',[.1 .1 .5 .8])


figure
scatter(foldcheck2(:,2),foldcheck2(:,1),6,'b','filled')
xlabel('Fold change, old code','fontsize',14)
ylabel('Fold change, new code','fontsize',14)
title('Fold Change, Perseus_enrichment_Gaussian_level_file.csv')
grid on


%% Make a figure or two explaining why example proteins are different in the new code


