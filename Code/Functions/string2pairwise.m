function string2pairwise(user)

%CORUM2PAIRWISE Creates a csv list of pairwise interactions
%    using the CORUM file allComplexes.csv.
%
%    Original script by Anders Kristensen.
%    Modified by Nichollas Scott, 2014.
%    Modified by Greg Stacey, 2016.



%% Read raw STRING file
fid = fopen(user.stringfile,'r');

% check header is correctly formatted
header = fgetl(fid);
if isempty(strfind(header,'protein1')) || isempty(strfind(header,'protein2'))
  error('corum2pairwise: Improperly formatted corum file (missing header?).')
end

% read body of file
cc = 0;
kk = 0;
pairwiselist = cell(10^6,2);
while 1
  t = fgetl(fid);
  if(~ischar(t)),break,end
  t1 = strsplit(t, ',');
  
  score = str2num(t1{4});
  if score>950
    cc = cc+1;
    pairwiselist{cc,1} = t1{1};
    pairwiselist{cc,2} = t1{2};
  else
    kk = kk+1;
  end
  
  if kk>100; break;  end
  
end
pairwiselist = pairwiselist(1:cc,:);
fclose(fid);



%% Write pairwise interactions to file

if ~exist([user.maindir '/Output'], 'dir'); mkdir([user.maindir '/Output']); end
if ~exist([user.maindir '/Output/tmp'], 'dir'); mkdir([user.maindir '/Output/tmp']); end
fn = [user.maindir '/Output/tmp/Corum_pairwise.csv'];
fid = fopen(fn,'w');
for ii = 1:size(pairwiselist,1)
  fprintf(fid,'%s,%s\n',pairwiselist{ii,1},pairwiselist{ii,2});
end
fclose(fid);


