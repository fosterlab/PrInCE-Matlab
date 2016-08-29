function out = readchromatogramfile2(fn)

% Makes:
%   num_val: first column replicate number, rest are fractions
%   txt_val: single column of protein IDs


fid = fopen(fn);

header = fgetl(fid);
header = strsplit(header,',');

cc = 0;
txt_val = cell(10^6,1);
txt_val{1} = header{1};
while ~feof(fid)
  t = fgetl(fid);
  t1 = strsplit(t,',');
  if cc ==0
    num_val = nan(10^6,length(t1)-1);
  end
  cc = cc+1;
  
  txt_val{cc+1} = t1{1};
  for kk = 1:size(num_val,2)
    num_val(cc,kk) = str2double(t1{kk+1});
  end
end

txt_val = txt_val(1:cc+1);
num_val = num_val(1:cc,:);

out.data = num_val;
out.textdata = txt_val;
