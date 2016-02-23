
%% Distribution of # of Gaussians

Ngauss = nan(Nchannels,Nproteins);
for ci = 1:Nchannels
  for ri = 1:Nproteins
    Ngauss(ci,ri) = length(Coef{ci,ri})/3;
  end
end

figure,hold on
h1 = hist(Ngauss(:),1:5);
bar(1:5,h1,'barwidth',0.9)
xlabel('Number of Gaussians used to model the chromatogram','fontsize',12)
ylabel('Count, number of chromatograms')
ylim([0 max(h1(:))*1.05])

set(gca,'xtick',1:5,'fontsize',12)

set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
sf=[figdir 'GaussBuild/Hist_NumberOfGaussians'];
saveas(gcf, sf, 'png');



%% Distribution of R^2

x = linspace(0,1,Nproteins/10);

figure,hold on
h1 = hist(adjrsquare(:),x);
bar(x,h1,'barwidth',0.9)
xlabel('R^2 between Guassian model and cleaned raw data','fontsize',12)
ylabel('Count, number of chromatograms')
xlim([-.01 1.01])

ax = axis;
x = ax(1) + diff(ax(1:2))*.35;
y1 = ax(3) + diff(ax(3:4))*.9;
y2 = ax(3) + diff(ax(3:4))*.8;
text(x,y1,['Average R^2 = ' num2str(nanmean(adjrsquare(:)))],'fontsize',12)
text(x,y2,['Median R^2 = ' num2str(nanmedian(adjrsquare(:)))],'fontsize',12)

set(gca,'fontsize',12)

set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
sf=[figdir 'GaussBuild/Hist_R2'];
saveas(gcf, sf, 'png');


%%

close all
