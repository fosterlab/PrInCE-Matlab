%MAKEFIGURES_GAUSSBUILD Makes figures for the PRINCE GaussBuild module.

n1 = size(rawdata{1},2);
n2 = size(cleandata{1},2);
frac1 = 1 + (n2-n1)/2;
frac2 = n2 - (n2-n1)/2;

xraw = 1:n1;
xclean = (1:n2) - (n2-n1)/2;


%% Heat map of raw chromatograms
try
  for ii = 1:Nchannels
    
    I = find(diff(replicate)>0);
    I = [I; length(replicate)];
    clear I2
    A = cell(length(I),1);
    A{1} = rawdata{ii}(1:I(1),:);
    [~,II] = max(A{1},[],2);
    [~,II2] = sort(II,'ascend');
    A{1} = A{1}(II2,:);
    ia = sum(not(isnan(A{1})),2)==0;
    A{1}(ia,:) = [];
    B = A{1};
    I2(1) = size(B,1);
    for jj = 2:length(I)
      A{jj} = rawdata{ii}(I(jj-1)+1:I(jj),:);
      [~,II] = max(A{jj},[],2);
      [~,II2] = sort(II,'ascend');
      A{jj} = A{jj}(II2,:);
      ia = sum(not(isnan(A{jj})),2)==0;
      A{jj}(ia,:) = [];
      B = [B; A{jj}];
      I2(jj) = size(B,1);
    end
    I = I2(1:end-1);
    
    figure
    imagesc(log10(B))
    axis xy
    hold on
    colormap bone
    
    x = xlim;
    y = ylim;
    text(x(1) + diff(x)*.8, y(1) + diff(y)*.05,'Replicate 1','color','w','fontsize',10)
    for jj = 1:length(I)
      plot(x, [1 1].*I(jj),'r')
      text(x(1) + diff(x)*.8, I(jj) + diff(y)*.05,['Replicate ' num2str(jj+1)],'color','w','fontsize',10)
    end
    
    xlabel('Fraction','fontsize',10)
    ylabel('Protein number','fontsize',10)
    title([user.silacratios{ii} ' - Raw'])
    
    set(gcf,'paperunits','inches','paperposition',[.25 2.5 6 9])
    sf=[figdir 'Chromatograms_raw_' user.silacratios{ii}];
    saveas(gcf, sf, 'png');
    
    ax = axis;
    
    I = find(diff(replicate)>0);
    I = [I; length(replicate)];
    clear I2
    A = cell(length(I),1);
    A{1} = cleandata{ii}(1:I(1),:);
    [~,II] = max(A{1},[],2);
    [~,II2] = sort(II,'ascend');
    A{1} = A{1}(II2,:);
    ia = sum(A{1}>0,2)==0;
    A{1}(ia,:) = [];
    B = A{1};
    I2(1) = size(B,1);
    for jj = 2:length(I)
      A{jj} = cleandata{ii}(I(jj-1)+1:I(jj),:);
      [~,II] = max(A{jj},[],2);
      [~,II2] = sort(II,'ascend');
      A{jj} = A{jj}(II2,:);
      ia = sum(A{jj}>0,2)==0;
      A{jj}(ia,:) = [];
      B = [B; A{jj}];
      I2(jj) = size(B,1);
    end
    I = I2(1:end-1);
    
    figure
    imagesc(log10(B(:,frac1:frac2)))
    axis xy
    hold on
    colormap bone
    
    x = xlim;
    y = ylim;
    text(x(1) + diff(x)*.8, y(1) + diff(y)*.05,'Replicate 1','color','w','fontsize',10)
    for jj = 1:length(I)
      plot(x, [1 1].*I(jj),'r')
      text(x(1) + diff(x)*.8, I(jj) + diff(y)*.05,['Replicate ' num2str(jj+1)],'color','w','fontsize',10)
    end
    
    xlabel('Fraction','fontsize',10)
    ylabel('Protein number','fontsize',10)
    title([user.silacratios{ii} ' - Clean'])
    
    set(gcf,'paperunits','inches','paperposition',[.25 2.5 6 9])
    sf=[figdir '/Chromatograms_clean_' user.silacratios{ii}];
    saveas(gcf, sf, 'png');
    
    axis(ax)
    
  end
catch
  warning('Failed to make heatmap figures.')
end


%% Distribution of # of Gaussians

Ngauss = nan(Nchannels,Nproteins);
for ci = 1:Nchannels
  for ri = 1:Nproteins
    Ngauss(ci,ri) = length(Coef{ci,ri})/3;
  end
end
Ngauss(Ngauss<1) = nan;

figure,hold on
h1 = hist(Ngauss(:),1:5);
bar(1:5,h1,'barwidth',0.9)
xlabel('Number of Gaussians used to model the chromatogram','fontsize',12)
ylabel('Count, number of chromatograms')
ylim([0 max(h1(:))*1.05])

set(gca,'xtick',1:5,'fontsize',12)

set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
sf=[figdir '/Hist_NumberOfGaussians'];
saveas(gcf, sf, 'png');



%% Distribution of Gaussian coefficients

Hall = cell(Nchannels,1);
Call = cell(Nchannels,1);
Wall = cell(Nchannels,1);
for ci = 1:Nchannels
  Hall{ci} = [];
  Call{ci} = [];
  Wall{ci} = [];
  for ri = 1:Nproteins
    H = Coef{ci,ri}(1:3:end);
    Hall{ci} = [Hall{ci} H];
    
    C = Coef{ci,ri}(2:3:end);
    Call{ci} = [Call{ci} C];
    
    W = Coef{ci,ri}(3:3:end);
    Wall{ci} = [Wall{ci} W];
  end
  
  % clean up parameters
  Hall{ci} = Hall{ci}(Hall{ci}>0 & Hall{ci}<250);
  Call{ci} = Call{ci}(Call{ci}>-10);
  Wall{ci} = Wall{ci}(Wall{ci}>0 & Wall{ci}<100);
end


figure,hold on
for ci = 1:Nchannels
  subplot(Nchannels,3,1 + (ci-1)*3)
  hist(Hall{ci},25)
  title(['H, ' user.silacratios{ci}])
  ylabel('Count','fontsize',10)
  xlabel('H, Height','fontsize',10)
  
  subplot(Nchannels,3,2 + (ci-1)*3)
  hist(Call{ci},25)
  title(['C, ' user.silacratios{ci}])
  ylabel('Count','fontsize',10)
  xlabel('C, Center fraction','fontsize',10)
  
  subplot(Nchannels,3,3 + (ci-1)*3)
  hist(Wall{ci},25)
  title(['W, ' user.silacratios{ci}])
  ylabel('Count','fontsize',10)
  xlabel('W, Width','fontsize',10)
end

set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
sf=[figdir '/Hist_GaussianParameters'];
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
sf=[figdir '/Hist_R2'];
saveas(gcf, sf, 'png');


%% Raw chromatograms, cleaned chromatograms, fits

if ~isfield(user,'fastgaussbuild')
  user.fastgaussbuild = 1;
end

if user.fastgaussbuild==1
  
  figdir1 = [figdir '/Chromatograms/'];
  if ~exist(figdir1, 'dir'); mkdir(figdir1); end
  
  colour_to_use=[0.254 0.411 0.882 %Colour: Royal Blue
    205/255 92/255 92/255 %Colour: Indian Red
    106/255 90/255 205/255 %Colour: Slate Blue
    135/255 206/255 250/255 %Colour: Sky Blue
    0.28 0.82 0.8 %Colour: Medium Turquoise
    178/255 34/255 34/255 %Colour: Firebrick
    255/255 69/255 0 %Colour: Orange Red
    244/255 238/255 224/255 %Colour: Honeydew 2
    34/255 139/255 34/255 %Colour: Forest Green
    222/255 184/255 135/255 %Colour: Burlywood
    186/255 85/255 211/255 %Colour: Medium Orchid
    219/255 112/255 147/255]; %Colour: Pale Violet Red
  
  h = figure;
  
  for ii = 1:Nproteins
    
    if sum(isnan(rawdata{1}(ii,:)))==length(rawdata{1}(ii,:))
      continue;
    end
    
    set(0, 'CurrentFigure', h);
    clf reset;
    
    % plot the raw data
    subplot(3,1,1),hold on
    for jj = 1:Nchannels
      plot(xraw,rawdata{jj}(ii,:),'s',...
        'MarkerFaceColor',colour_to_use(jj,:),'MarkerSize',4);
    end
    for jj = 1:Nchannels
      plot(xraw,rawdata{jj}(ii,:),'color',colour_to_use(jj,:));
    end
    xlim([xraw(1)-1 xraw(end)+1])
    y = ylim;
    if y(1)>0
      y(1) = 0.01;
    end
    ylim(y)
    xlabel('Fraction','fontsize',10)
    ylabel('Isotopologue ratio','FontSize', 10);
    legend(user.silacratios,'location','best')
    title('Raw chromatogram','fontsize',10)
    ax = axis;
    
    % plot the clean data
    subplot(3,1,2),hold on
    for jj = 1:Nchannels
      plot(xclean,cleandata{jj}(ii,:),'s',...
        'MarkerFaceColor',colour_to_use(jj,:),'MarkerSize',4);
      plot(xclean,cleandata{jj}(ii,:),'color',colour_to_use(jj,:));
    end
    xlabel('Fraction','fontsize',10)
    ylabel('Isotopologue ratio','FontSize', 10);
    title('Cleaned chromatogram','fontsize',10)
    axis(ax)
    
    % plot the Gaussians
    subplot(3,1,3),hold on
    xfit = linspace(xraw(1)-1,xraw(end)+1,101);
    yfit = zeros(Nchannels,length(xfit));
    Ngauss = nan(Nchannels,1);
    for jj = 1:Nchannels
      Ngauss(jj) = length(Coef{jj,ii})/3;
      if Ngauss(jj) == 0; continue;end
      for kk = 1:Ngauss(jj)
        H = Coef{jj,ii}(1 + (kk-1)*3);
        C = Coef{jj,ii}(2 + (kk-1)*3);
        W = Coef{jj,ii}(3 + (kk-1)*3);
        yfit(jj,:) = yfit(jj,:) + H*exp(-((xfit-C)/W).^2);
      end
      plot(xfit,yfit(jj,:),'color',colour_to_use(jj,:),'linewidth',2);
      if Ngauss(jj)>0
        s1 = [num2str(Ngauss(jj)) ' Gauss, R^2=' num2str(round(adjrsquare(jj,ii)*100)/100)];
        text(ax(1)+diff(ax(1:2))*.75, ax(3)+diff(ax(3:4))*(1-.1*jj), s1)
      end
    end
    if sum(Ngauss(:))>0
      legend(user.silacratios,'location','northwest')
    end
    xlabel('Fraction','fontsize',10)
    ylabel('Isotopologue ratio','FontSize', 10);
    title('Fitted Gaussians','fontsize',10)
    axis(ax)
    
    
    sf = [figdir '/Chromatograms/' mat2str(replicate(ii)),'_',txt_val{1}{ii+1,1},'.png'];
    saveas(gcf, sf);
    
  end
end

%%

close all
