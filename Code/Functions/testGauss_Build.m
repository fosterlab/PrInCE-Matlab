% Idea: Fit all the chromatograms with an AIC/AICc/BIC method and compare to Nick's fitting.


% Load in Nick's Gaussian fits
fn = '/Users/Mercy/Academics/Foster/NickCodeData/1_Gaussian processing/MvsL/MvsL_Combined_OutputGaus.csv';
fn2 = fopen(fn);
%Processed_data = textscan(fn2, '%n %n %n %s %n %n %n %n %n %n', 'Delimiter', ',');
data = importdata(fn);

Protein_number = zeros(size(data.data,1)-1,1);
replicate = zeros(size(data.data,1)-1,1);
Protein_name = cell(size(data.data,1)-1,1);
for ii=2:size(data.data,1)
  Protein_number(ii-1) = str2num(data.textdata{ii,2});
  replicate(ii-1) = str2num(data.textdata{ii,3});
  Protein_name{ii-1} = data.textdata{ii,4};
end

Nprot = max(Protein_number);
goodprot = unique(Protein_number);

f = @(x,a,b,c) a(1) * exp(-(x - b(1)).^2 / c(1)^2) + a(2) * exp(-(x - b(2)).^2 / c(2)^2) + ...
  a(3) * exp(-(x - b(3)).^2 / c(3)^2) + a(4) * exp(-(x - b(4)).^2 / c(4)^2) + ...
  a(5) * exp(-(x - b(5)).^2  / c(5)^2);
x = (1:65)' - 5;

cols = [1 1 0; 1 0 1; 0 1 1; .6 .6 .6; 1 .5 .5];

for ri = goodprot(500:end)'
  rawchrom = rawdata{1}(ri,:);
  cleanchrom = cleandata{1}(ri,:);
  
  % Nick's fit
  a = zeros(5,1);
  b = zeros(5,1);
  c = zeros(5,1);
  I = find(Protein_number == ri);
  for gi = 1:length(I)
    a(gi) = data.data(I(gi),1);
    b(gi) = data.data(I(gi),2);
    c(gi) = data.data(I(gi),3);
  end
  yhat_nick = f(x,a,b,c);
  
  % My fit
  a2 = zeros(5,1);
  b2 = zeros(5,1);
  c2 = zeros(5,1);
  model = choosemodel_AIC(cleandata{1}(ri,:),x,'AICc');
  for gi = 1:length(model.coeffs)/3
    a2(gi) = model.coeffs((gi-1)*3+1);
    b2(gi) = model.coeffs((gi-1)*3+2);
    c2(gi) = model.coeffs((gi-1)*3+3);
  end
  yhat_greg = f(x,a2,b2,c2);
  
  figure
  subplot(2,1,1),hold on
  plot(rawchrom,'k')
  plot(x,cleanchrom,'color',[.6 .6 .6],'linewidth',2)
  plot(x,yhat_nick,'r')
  for gi = 1:length(I)
    plot(x',f(x,[a(gi) 0 0 0 0],[b(gi) 0 0 0 0],[c(gi) 0 0 0 0]),'color',cols(gi,:),'linestyle','--')
  end
  
  subplot(2,1,2),hold on
  plot(rawchrom,'k')
  plot(x,cleanchrom,'color',[.6 .6 .6],'linewidth',2)
  plot(x,yhat_greg,'g')
  for gi = 1:length(model.coeffs)/3
    plot(x',f(x,[a2(gi) 0 0 0 0],[b2(gi) 0 0 0 0],[c2(gi) 0 0 0 0]),'color',cols(gi,:),'linestyle','--')
  end
  pause
end
