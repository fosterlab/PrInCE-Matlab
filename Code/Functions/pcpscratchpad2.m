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



%%

% noise, 1000 random events between 0 and 10,000 ms
licks = rand(1,500) * 10000;
% signal, events evenly spaced with a period of 1/(60 Hz)
electrical_noise = 0 : 1/60 * 1000 : 10000;
electrical_noise(randsample(1:length(electrical_noise),round(length(electrical_noise)*.10))) = [];
% Simulated data = licks + electrical_noise
X = sort([licks electrical_noise]);

% inter-lick intervals
ILI = diff(X);

figure
x = linspace(0,max(X),51);
h = hist(X,x);
bar(x,h,.8)
xlim([min(X) max(X)])
xlabel('Time (ms)')
ylabel('Count')
title('Histogram of licks + electrical noise')
set(gcf,'units','normalized','position',[.1 .1 .7 .4])

figure,hold on
x = 0 : 1/600*1000 : max(ILI);
h = hist(ILI,x);
bar(x,h,.8)
y = ylim;
plot([1 1]*1/60*1000, y, '--r')
xlim([min(ILI) max(ILI)])
xlabel('Inter-lick interval (ms)')
ylabel('Count')
title('Histogram of ILIs')
set(gcf,'units','normalized','position',[.1 .1 .7 .4])

% dbin_lick=10; 
% binrange=0:dbin_lick8max(x); 
% bin_count=histc(X,binrange); 
% bar(binrange,bin_count,'histc'); 
% Fs = 1000/dbin_lick; % Sampling frequency 
% T = X(end); % Sampling period 
% L = length(X); % Length of signal 
% t = (0:L-1)*T; 
% Y = fft(X); 
% P2 = abs(Y/L); 
% P1 = P2(1:L/2+1); 
% P1(2:end-1) = 2*P1(2:end-1); 
% %f = Fs*(0:(L/2))/L;
% f = (0:(L/2)) * Fs / L;
% figure
% plot(f,P1) 
% title('Single-Sided Amplitude Spectrum of X(t)') 
% 
% 
% xlabel('f (Hz)') 
% ylabel('|P1(f)|')

