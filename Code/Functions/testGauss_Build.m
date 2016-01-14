% Idea: Fit all the chromatograms with an AIC/AICc/BIC method and compare to Nick's fitting.
%
% Figures to make:
% Fig. 1a: Venn diagram of proteins fit with at least one Gaussian by i) me and ii) Nick. How much 
%   overlap is there? Does one of us just fit more Gaussians?
% Fig. 1b: For each protein with a Gaussian fit, reconstruct the chromatogram from the Gaussian 
%   parameters. Then calculate a correlation coefficient between mine and Nick?s data. How similarly 
%   do we model each chromatogram?

% Load in Nick's Gaussian fits
fn = '/Users/Mercy/Academics/Foster/NickCodeData/1_Gaussian processing/MvsL/MvsL_Combined_OutputGaus.csv';
data = importdata(fn);
Nprot = size(data.data,1);
Protein_number = zeros(size(data.data,1)-1,1);
replicate = zeros(size(data.data,1)-1,1);
Protein_name = cell(size(data.data,1)-1,1);
Height = zeros(Nprot,1);
Center = zeros(Nprot,1);
Width = zeros(Nprot,1);
for ii = 2:Nprot+1
  Protein_number(ii-1) = str2num(data.textdata{ii,2});
  replicate(ii-1) = str2num(data.textdata{ii,3});
  Protein_name{ii-1} = data.textdata{ii,4};
  Height(ii-1) = data.data(ii-1,1);
  Center(ii-1) = data.data(ii-1,2);
  Width(ii-1) = data.data(ii-1,3);
end

% Load in my Gaussian fits
fn = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/Data/GaussBuild/MvsL_Combined_OutputGaus.csv';
datag = importdata(fn);
Nprot = size(datag.data,1);
Protein_numberg = zeros(size(datag.data,1)-1,1);
replicateg = zeros(size(datag.data,1)-1,1);
Protein_nameg = cell(size(datag.data,1)-1,1);
Heightg = zeros(Nprot,1);
Centerg = zeros(Nprot,1);
Widthg = zeros(Nprot,1);
for ii = 2:Nprot+1
  Protein_numberg(ii-1) = str2num(datag.textdata{ii,2});
  replicateg(ii-1) = str2num(datag.textdata{ii,3});
  Protein_nameg{ii-1} = datag.textdata{ii,4};
  Heightg(ii-1) = datag.data(ii-1,1);
  Centerg(ii-1) = datag.data(ii-1,2);
  Widthg(ii-1) = datag.data(ii-1,3);
end

f = @(x,a,b,c) a(1) * exp(-(x - b(1)).^2 / c(1)^2) + a(2) * exp(-(x - b(2)).^2 / c(2)^2) + ...
  a(3) * exp(-(x - b(3)).^2 / c(3)^2) + a(4) * exp(-(x - b(4)).^2 / c(4)^2) + ...
  a(5) * exp(-(x - b(5)).^2  / c(5)^2);
x = (1:65)' - 5;
x(5) = 0.000000001;

cols = [1 1 0; 1 0 1; 0 1 1; .6 .6 .6; 1 .5 .5]*.8;

% for each of Nick's protein
%   get all of Nick's Gaussian parameters (a,b,c)
%   get all of my Gaussian parameters (a2,b2,c2)
%   get the raw chromatogram
goodprot = unique(Protein_number);
Ng_nick = zeros(length(goodprot),1);
Ng_greg = zeros(length(goodprot),1);
R2_ng = zeros(length(goodprot),1);
yhat_nick = zeros(length(goodprot),65);
yhat_greg = zeros(length(goodprot),65);
cleanchrom = zeros(length(goodprot),65);
sumA = zeros(length(goodprot),1);
for ri = 1:length(goodprot)
  
  % Nick's fit
  In = find(Protein_number==goodprot(ri));
  a = zeros(5,1);
  b = zeros(5,1);
  c = zeros(5,1);
  for gi = 1:length(In)
    a(gi) = Height(In(gi));
    b(gi) = Center(In(gi));
    c(gi) = Width(In(gi));
  end
  yhat_nick(ri,:) = f(x,a,b,c);
  Ng_nick(ri) = length(In);
  
  % My fit
  Ig = find(Protein_numberg==goodprot(ri));
  a2 = zeros(5,1);
  b2 = zeros(5,1);
  c2 = zeros(5,1);
  for gi = 1:length(Ig)
    a2(gi) = Heightg(Ig(gi));
    b2(gi) = Centerg(Ig(gi));
    c2(gi) = Widthg(Ig(gi));
  end
  yhat_greg(ri,:) = f(x,a2,b2,c2);
  Ng_greg(ri) = length(Ig);
  
  tmp = corrcoef(yhat_nick(ri,:),yhat_greg(ri,:));
  R2_ng(ri) = tmp(1,2)^2;
  %if R2_ng(ri)<0.1;a2,pause;end
  sumA(ri) = sum(a2);
  
  % Raw chromatogram
  pn = Protein_numberg(Ig);
  cleanchrom(ri,:) = cleandata{1}(pn(1),:);
  
%   % Raw chromatogram
%   pn = Protein_numberg(Ig);
%   rawchrom = rawdata{1}(pn,:);
%   
%   figure
%   subplot(2,1,1),hold on
%   plot(rawchrom','k')
%   plot(x,yhat_nick,'r')
%   for gi = 1:length(In)
%     plot(x',f(x,[a(gi) 0 0 0 0],[b(gi) 0 0 0 0],[c(gi) 0 0 0 0]),'color',cols(gi,:),'linestyle','--')
%   end
%   
%   subplot(2,1,2),hold on
%   plot(rawchrom','k')
%   plot(x,yhat_greg,'g')
%   for gi = 1:length(Ig)
%     plot(x',f(x,[a2(gi) 0 0 0 0],[b2(gi) 0 0 0 0],[c2(gi) 0 0 0 0]),'color',cols(gi,:),'linestyle','--')
%   end
%   set(gcf,'units','normalized','position',[.1 .1 .3 .8])
%   
%   pause
end




%% Venn diagram

Allprots = unique([unique(Protein_name);unique(Protein_nameg)]);
Protg = find(ismember(Allprots,Protein_nameg));
Protn = find(ismember(Allprots,Protein_name));
figure
myVenn2([length(Protn) length(Protg)], length(Allprots))
%text(0,0,num2str(length(Inc2)))
%text(-15,10,num2str(length(Incn) - length(Inc2)))
%text(16,10,num2str(length(Incg) - length(Inc2)))
set(gca,'xtick',[],'ytick',[])




%% Number-of-Gaussians distribution

figure
subplot(2,1,1),hold on
h1 = hist(Ng_nick,1:5);
bar(h1,'barwidth',0.9)
text(4,2000,'Old code (holdout)')
xlabel('Number of Gaussians fit, old code','fontsize',12)
ylabel('Count')

subplot(2,1,2),hold on
h2 = hist(Ng_greg,1:5);
bar(h2,'barwidth',0.9)
text(4,2100,'New code (robust + AIC)')
xlabel('Number of Gaussians fit, new code','fontsize',12)
ylabel('Count')



%% Examples of good, okay, bad

goodI = find(R2_ng>.95);
goodI = randsample(goodI,1);
okayI = find(R2_ng<.7 & R2_ng>0.4);
okayI = randsample(okayI,1);
badI = find(R2_ng<.2);
badI = randsample(badI,1);
%goodI = 4671;
%okayI = 1306;
%badI = 761;

figure
subplot(2,2,1)
hist(R2_ng,linspace(0,1,101))
xlim([-.01 1.01])
xlabel('R^2 between old and new Gaussian curves','fontsize',12)
ylabel('Count','fontsize',12)
subplot(2,2,2),hold on
plot(cleanchrom(goodI,:),'color',[.7 .7 .7],'linewidth',2)
plot(yhat_nick(goodI,:),'k')
plot(yhat_greg(goodI,:),'g')
legend('Chromat.','Old Gauss','New Gauss','location','best')
y = ylim;
text(35,y(1) + diff(y)*.75,['R^2 = ' num2str(R2_ng(goodI))])
xlim([0 65])
xlabel('Fraction','fontsize',12)
subplot(2,2,3),hold on
plot(cleanchrom(okayI,:),'color',[.7 .7 .7],'linewidth',2)
plot(yhat_nick(okayI,:),'k')
plot(yhat_greg(okayI,:),'g')
y = ylim;
text(35,y(1) + diff(y)*.75,['R^2 = ' num2str(R2_ng(okayI))])
xlim([0 65])
xlabel('Fraction','fontsize',12)
subplot(2,2,4),hold on
plot(cleanchrom(badI,:),'color',[.7 .7 .7],'linewidth',2)
plot(yhat_nick(badI,:),'k')
plot(yhat_greg(badI,:),'g')
y = ylim;
text(35,y(1) + diff(y)*.75,['R^2 = ' num2str(R2_ng(badI))])
xlim([0 65])
xlabel('Fraction','fontsize',12)

set(gcf,'units','normalized','position',[.1 .1 .5 .7])


