% Compares the output of ROC_PCPSILAC.m to Nick's old code. Specifically, does the new code find the
% same interactions? Treating Nick's interactions as TP, what is the new code's precision/recall?

spl = 0;

% % Load in Nick's interactions
% fn = '/Users/Mercy/Academics/Foster/NickCodeData/4B_ROC_homologue_DB/Combined Cyto new DB/Combined results/Interactions_across_replicate_70pc.csv';
% dataOld = importdata(fn);
% dataOld.textdata(1) = [];
%
% % Load in Nick's interactions
% fn = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/Data/ROC/CombinedResults/Interactions_across_replicate_70pcb.csv';
% dataNew = importdata(fn);
% dataNew.textdata(1) = [];

% Load in New interactions
%fn = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/Data/ROC/CombinedResults/Final_Interactions_list_70_precisionb.csv';
%fn = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/Data/ROC/CombinedResults/Final_Interactions_list_47_precisionb.csv';
%fn = '/Users/Mercy/Academics/Foster/Tissue_PCPSILAC/PCPSILAC_analysis/Data/ROC/CombinedResults/Final_Interactions_list_47_precision.csv';
%fn = '/Users/Mercy/Academics/Foster/NickCodeData/Old runs/GregPCP_20160517/Data/ROC/CombinedResults/Final_Interactions_list_70_precision.csv';
fn = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/Data/Output/Interactions/Final_Interactions_list_47_precision.csv';
dataNew.data = zeros(20000,9);
dataNew.text = cell(20000,3);
j = 0;
fid = fopen(fn,'r');
fgetl(fid); % skip header
while 1
  t = fgetl(fid);
  if(~ischar(t)),break,end
  j = j+1;
  t1 = strsplit(t, ',');
  dataNew.data(j,1) = str2num(t1{10}); % proteins in corum?
  dataNew.data(j,2) = str2num(t1{11}); % interaction in corum?
  for ii=1:6
    dataNew.data(j,2+ii) = ismember(num2str(ii),t1{6}); % interaction in corum?
  end
  dataNew.data(j,9) = str2num(t1{13});
  dataNew.text{j,1} = t1{1};
  dataNew.text{j,2} = t1{2};
  dataNew.text{j,3} = t1{3};
end
dataNew.data = dataNew.data(1:j,:);
dataNew.text = dataNew.text(1:j,:);
% Load in Old interactions
fn = '/Users/Mercy/Academics/Foster/NickCodeData/4B_ROC_homologue_DB/Combined Cyto new DB/Combined results/Final_Interactions_list_70_precision.csv';
%fn = '/Users/Mercy/Downloads/higer precision tissue data/Combined results/Final_Interactions_list_80_precision.csv';
dataOld.data = zeros(20000,8);
dataOld.text = cell(20000,3);
j = 0;
fid = fopen(fn,'r');
fgetl(fid); % skip header
while 1
  t = fgetl(fid);
  if(~ischar(t)),break,end
  j = j+1;
  t1 = strsplit(t, ',');
  dataOld.data(j,1) = str2double(t1{11}); % proteins in corum?
  dataOld.data(j,2) = str2double(t1{12}); % interaction in corum?
  for ii=1:6
    dataOld.data(j,2+ii) = ismember(num2str(ii),t1{6}); % interaction in corum?
  end
  dataOld.text{j,1} = t1{1};
  dataOld.text{j,2} = t1{2};
  dataOld.text{j,3} = t1{3};
end
dataOld.data = dataOld.data(1:j,:);
dataOld.text = dataOld.text(1:j,:);


%% Remove duplicate interactions (just keep protein-protein pairs)

dupOld = zeros(size(dataOld.text,1),1);
for ii = 1:size(dataOld.text,1)
  dupOld(ii) = sum(strcmp(dataOld.text(ii,1),dataOld.text(:,1)));
end
dupNew = zeros(size(dataNew.text,1),1);
for ii = 1:size(dataNew.text,1)
  dupNew(ii) = sum(strcmp(dataNew.text(ii,1),dataNew.text(:,1)));
end

dataOld.data = dataOld.data(dupOld==1,:);
dataOld.text = dataOld.text(dupOld==1,:);
dataNew.data = dataNew.data(dupNew==1,:);
dataNew.text = dataNew.text(dupNew==1,:);


%% For each replicate, what's the interaction overlap between old and new?
figure
for ii = 1:6
  
  I = dataOld.data(:,ii+2)>0; % interaction found in this replicate
  intold = dataOld.text(I,1);
  
  I = dataNew.data(:,ii+2)>0; % interaction found in this replicate
  intnew = dataNew.text(I,1);
  
  % For each interaction in Nick's code, see if it's in my code
  isin1 = zeros(size(intold));
  for jj = 1:length(intold)
    I = findstr(intold{jj},'_');
    protA = intold{jj}(1:I-1);
    protB = intold{jj}(I+1:end);
    
    isin1(jj) = sum(strcmp([protA '_' protB],intnew)) + sum(strcmp([protB '_' protA],intnew));
  end
  if sum(isin1(jj)>1)>0;1,pause;end
  
  % For each interaction in my code, see if it's in Nick's code
  isin2 = zeros(size(intnew));
  for jj = 1:length(intnew)
    I = findstr(intnew{jj},'_');
    protA = intnew{jj}(1:I-1);
    protB = intnew{jj}(I+1:end);
    
    isin2(jj) = sum(strcmp([protA '_' protB],intold)) + sum(strcmp([protB '_' protA],intold));
  end
  
  if sum(isin1>0)~=sum(isin2);disp('uh oh111');end
  
  nnew = length(isin2);
  nold = length(isin1);
  
  TP = sum(isin1);
  FP = sum(isin2==0);
  FN = sum(isin1==0);
  prec = round(TP/(TP+FP)*100)/100;
  rec = round(TP/(TP+FN)*100)/100;
  
  subplot(2,3,ii)
  myVenn2([nnew nold], TP)
  axis([-50 100 -55 55])
  ax = axis;
  text(ax(1)+0.05*diff(ax(1:2)),ax(3)+0.95*diff(ax(3:4)),['prec=' num2str(prec)])
  text(ax(1)+0.05*diff(ax(1:2)),ax(3)+0.9*diff(ax(3:4)),['rec=' num2str(rec)])
  text(ax(1)+0.6*diff(ax(1:2)),ax(3)+0.95*diff(ax(3:4)),['N=' num2str(length(isin2))],'color',[1 0 0])
  text(ax(1)+0.6*diff(ax(1:2)),ax(3)+0.9*diff(ax(3:4)),['N=' num2str(length(isin1))],'color',[0 1 0])
  title(['Interactions in replicate ' num2str(ii)])
  axis(ax)
  
  if ii==1
    legend('New','Old','location','southwest')
  end
end
pause(.001)

if spl
  set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
  graphdir = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/Figures/Test/';
  saveas(gcf,[graphdir 'ROC_fig1'],'jpg')
end


%%
% For interactions found in at least N replicates, what's the interaction overlap between old and new?
figure
for ii = 1:6
  
  I = sum(dataOld.data(:,3:8)>0,2)>=ii; % interaction found in at least ii replicates
  intold = dataOld.text(I,1);
  
  I = sum(dataNew.data(:,3:8)>0,2)>=ii; % interaction found in at least ii replicates
  intnew = dataNew.text(I,1);
  
  % For each interaction in Nick's code, see if it's in my code
  isin1 = zeros(size(intold));
  for jj = 1:length(intold)
    I = findstr(intold{jj},'_');
    protA = intold{jj}(1:I-1);
    protB = intold{jj}(I+1:end);
    
    isin1(jj) = sum(strcmp([protA '_' protB],intnew)) + sum(strcmp([protB '_' protA],intnew));
  end
  
  % For each interaction in my code, see if it's in Nick's code
  isin2 = zeros(size(intnew));
  for jj = 1:length(intnew)
    I = findstr(intnew{jj},'_');
    protA = intnew{jj}(1:I-1);
    protB = intnew{jj}(I+1:end);
    
    isin2(jj) = sum(strcmp([protA '_' protB],intold)) + sum(strcmp([protB '_' protA],intold));
  end
  
  if sum(isin1>0)~=sum(isin2);disp('uh oh111');end
  
  nnew = length(isin2);
  nold = length(isin1);
  
  TP = sum(isin1);
  FP = sum(isin2==0);
  FN = sum(isin1==0);
  prec = round(TP/(TP+FP)*100)/100;
  rec = round(TP/(TP+FN)*100)/100;
  
  subplot(2,3,ii)
  myVenn2([nnew nold], TP)
  axis square equal
  ax = axis;
  text(ax(1)+0.05*diff(ax(1:2)),ax(3)+0.95*diff(ax(3:4)),['prec=' num2str(prec)])
  text(ax(1)+0.05*diff(ax(1:2)),ax(3)+0.9*diff(ax(3:4)),['rec=' num2str(rec)])
  text(ax(1)+0.6*diff(ax(1:2)),ax(3)+0.95*diff(ax(3:4)),['N=' num2str(length(isin2))],'color',[1 0 0])
  text(ax(1)+0.6*diff(ax(1:2)),ax(3)+0.9*diff(ax(3:4)),['N=' num2str(length(isin1))],'color',[0 1 0])
  title(['Interactions in >=' num2str(ii) ' replicates'])
  set(gca,'xtick',-75:10:150,'ytick',-70:10:70,'xticklabel','','yticklabel','')
  axis([ax(1) ax(2) ax(3)-diff(ax(3:4))*.05 ax(4)+diff(ax(3:4))*.05])
  grid on
  pause(.001)
  if ii==1
    legend('New','Old','location','southwest')
  end
end
pause(.001)

if spl
  set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
  graphdir = '/Users/Mercy/Desktop/';
  saveas(gcf,[graphdir 'ROC_fig2'],'jpg')
end


%% Calculate final precision/recall

Precold = zeros(6,2);
Precnew = zeros(6,2);
for ii=1:6
  I = dataOld.data(:,1)==1 & sum(dataOld.data(:,3:8),2)>=ii;
  TP = sum(I & dataOld.data(:,2)>0);
  FP = sum(I & dataOld.data(:,2)==0);
  Precold(ii,1) = TP/(TP+FP);
  Precold(ii,2) = TP+FP;
  
  I = dataNew.data(:,1)==1 & sum(dataNew.data(:,3:8),2)>=ii;
  TP = sum(I & dataNew.data(:,2)==1);
  FP = sum(I & dataNew.data(:,2)==0);
  Precnew(ii,1) = TP/(TP+FP);
  Precnew(ii,2) = TP+FP;
end

figure
subplot(2,1,1),hold on
plot(Precold(:,1),'g')
plot(Precnew(:,1),'r')
legend('old','new')
set(gca,'xtick',1:6)
grid on
subplot(2,1,2),hold on
scatter(Precold(:,1),Precold(:,2),'g','filled')
scatter(Precnew(:,1),Precnew(:,2),'r','filled')

if spl
  set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
  graphdir = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/Figures/Test/';
  saveas(gcf,[graphdir 'ROC_fig3'],'jpg')
end


Precold = zeros(6,2);
Precnew = zeros(6,2);
for ii=1:6
  I = dataOld.data(:,1)==1 & dataOld.data(:,ii+2);
  TP = sum(I & dataOld.data(:,2)>0);
  FP = sum(I & dataOld.data(:,2)==0);
  Precold(ii,1) = TP/(TP+FP);
  Precold(ii,2) = TP+FP;
  
  I = dataNew.data(:,1)==1 & dataNew.data(:,ii+2);
  TP = sum(I & dataNew.data(:,2)==1);
  FP = sum(I & dataNew.data(:,2)==0);
  Precnew(ii,1) = TP/(TP+FP);
  Precnew(ii,2) = TP+FP;
end

figure
subplot(2,1,1),hold on
bar([Precold(:,1) Precnew(:,1)])
x = xlim;
plot(x,[0.7 0.7],'--r')
ylim([0 1])
set(gca,'xtick',1:6)
grid on
subplot(2,1,2),hold on
scatter(Precold(:,1),Precold(:,2),'g','filled')
scatter(Precnew(:,1),Precnew(:,2),'r','filled')
pause(.001)

if spl
  set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
  graphdir = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/Figures/Test/';
  saveas(gcf,[graphdir 'ROC_fig4'],'jpg')
end



%% Nick's bar graph of i) FP, ii) TP, iii) proteins not in CORUM

fpold = zeros(6,1);
tpold = zeros(6,1);
freshold = zeros(6,1);
fpnew = zeros(6,1);
tpnew = zeros(6,1);
freshnew = zeros(6,1);
for ii=1:6
  I = dataOld.data(:,1)==1;
  I2 = sum(dataOld.data(:,3:8),2)>=ii;
  tpold(ii) = sum(I & I2 & dataOld.data(:,2)>0);
  fpold(ii) = sum(I & I2 & dataOld.data(:,2)==0);
  freshold(ii) = sum(~I & I2);
  
  I = dataNew.data(:,1)==1;
  I2 = sum(dataNew.data(:,3:8),2)>=ii;
  tpnew(ii) = sum(I & I2 & dataNew.data(:,2)==1);
  fpnew(ii) = sum(I & I2 & dataNew.data(:,2)==0);
  freshnew(ii) = sum(~I & I2);
end

a = [fpold tpold freshold];
b = [fpnew tpnew freshnew];
x = 1:6;

figure
subplot(2,1,1),hold on
B1 = bar(x-0.125,a,0.25,'stacked');
B2 = bar(x+0.125,b,0.25,'stacked');
C1=[.5 .5 1; 1/sqrt(2) 1/sqrt(2) .5; 1 .5 .5];
C2=[0 0 1; 1/sqrt(2) 1/sqrt(2) 0; 1 0 0] * .7;
for n=1:length(B1)
  set(B1(n),'facecolor',C1(n,:));
end
for n=1:length(B2)
  set(B2(n),'facecolor',C2(n,:));
end
y = ylim;
text(1-.3,sum(a(1,:))+diff(y)*.03,'Old')
text(1+.05,sum(b(1,:))+diff(y)*.03,'New')
legend('FP','TP','New')
grid on
ylabel('Count')
title('Number of interactions found in at least N replicates')
set(gca,'xtick',1:6,'xticklabel',{'>=1' '>=2' '>=3' '>=4' '>=5' '6'})

subplot(2,1,2),hold on
B1 = bar(x-0.125,a,0.25,'stacked');
B2 = bar(x+0.125,b,0.25,'stacked');
C1=[.5 .5 1; 1/sqrt(2) 1/sqrt(2) .5; 1 .5 .5];
C2=[0 0 1; 1/sqrt(2) 1/sqrt(2) 0; 1 0 0] * .7;
for n=1:length(B1)
  set(B1(n),'facecolor',C1(n,:));
end
for n=1:length(B2)
  set(B2(n),'facecolor',C2(n,:));
end
ylim([0 ceil(b(2,2)/200)*200])
y = ylim;
text(6-.3,sum(a(6,:))+diff(y)*.05,'Old')
text(6+.05,sum(b(6,:))+diff(y)*.05,'New')
legend('FP','TP','New')
grid on
xlabel('Number of replicates')
ylabel('Count')
set(gca,'xtick',1:6,'xticklabel',{'>=1' '>=2' '>=3' '>=4' '>=5' '6'})

if spl
  set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
  graphdir = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/Figures/Test/';
  saveas(gcf,[graphdir 'ROC_fig5'],'jpg')
end


%% Calculate precision/recall for each replicate
% This is really just a check that they're all at the desired level
if 0
  for rep = 1:6
    sf = [datadir2 'Interaction_list_pc70_rep' num2str(rep) '.mat'];
    load(sf)
    
    clear FNBinaryInteractions List_of_final_interactions_Protein I I2
    for ii = 1:size(BinaryInteractions,1)
      I(ii) = BinaryInteractions{ii,6};
      I2(ii) = BinaryInteractions{ii,7};
    end
    TP = sum(I & I2);
    FP = sum(I & ~I2);
    prec(rep) = TP/(TP+FP);
    length(I)
  end
  figure,hold on
  bar(prec)
  x = xlim;
  plot(x,[0.7 0.7],'--r')
  xlim(x)
  ylim([0 1])
  
  if spl
    set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
    graphdir = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/Figures/Test/';
    saveas(gcf,[graphdir 'ROC_fig7'],'jpg')
  end
end


%% Visualize TP,FP,TN,FN

tp = find(Final_interactions & TP_Matrix & inverse_self);
fp = find(Final_interactions & ~TP_Matrix & inverse_self);
tn = find(~Final_interactions & ~TP_Matrix & inverse_self);
fn = find(~Final_interactions & TP_Matrix & inverse_self);
a = randsample(length(tp),length(tp));
tp = tp(a);
[tp1, tp2] = ind2sub(size(TP_Matrix),tp);

a = randsample(length(fp),length(fp));
fp = fp(a);
[fp1, fp2] = ind2sub(size(TP_Matrix),fp);

a = randsample(length(tn),length(tn));
tn = tn(a);
[tn1, tn2] = ind2sub(size(TP_Matrix),tn);

a = randsample(length(fn),length(fn));
fn = fn(a);
[fn1, fn2] = ind2sub(size(TP_Matrix),fn);

figure
subplot(2,2,1),hold on
plot(Chromatograms(tp1(1),:),'r')
plot(Chromatograms(tp2(1),:),'g')
title('TP')
subplot(2,2,2),hold on
plot(Chromatograms(fp1(1),:),'r')
plot(Chromatograms(fp2(1),:),'g')
title('FP')
subplot(2,2,3),hold on
plot(Chromatograms(fn1(1),:),'r')
plot(Chromatograms(fn2(1),:),'g')
title('FN')
subplot(2,2,4),hold on
plot(Chromatograms(tn1(1),:),'r')
plot(Chromatograms(tn2(1),:),'g')
title('TN')


figure
for ii = 1:8
  subplot(4,8,ii),hold on
  plot(Chromatograms(tp1(ii),:),'r')
  plot(Chromatograms(tp2(ii),:),'g')
  if ii==1;ylabel('TP');end
  
  subplot(4,8,ii+8),hold on
  plot(Chromatograms(fp1(ii),:),'r')
  plot(Chromatograms(fp2(ii),:),'g')
  if ii==1;ylabel('FP');end
  
  subplot(4,8,ii+8*2),hold on
  plot(Chromatograms(fn1(ii),:),'r')
  plot(Chromatograms(fn2(ii),:),'g')
  if ii==1;ylabel('FN');end
  
  subplot(4,8,ii+8*3),hold on
  plot(Chromatograms(tn1(ii),:),'r')
  plot(Chromatograms(tn2(ii),:),'g')
  if ii==1;ylabel('TN');end
end
set(gcf,'units','normalized','position',[.05 .025 .9 .95])
set(gcf,'paperunits','normalized','paperposition',[.1 .1 .8 .5])



%% Visually check a few interactions

replicate_counter = 5;
sf = [datadir2 'score_rep' num2str(replicate_counter) '.mat'];
load(sf)
scoreMatrix = reshape(scoreMatrix,length(Protein.Center),length(Protein.Center));

prot1i = find(strcmp(Protein.Isoform,'O43660-2'));
prot2i = find(strcmp(Protein.Isoform,'P36578'));

figure,hold on
plot(Chromatograms(prot1i(1),:))
plot(Chromatograms(prot2i(1),:))
title(['Score = ' num2str(scoreMatrix(prot1i(1),prot2i(1)))])
xlim([0 65])


