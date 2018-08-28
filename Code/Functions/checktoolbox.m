function checktoolbox()

% CHECKTOOLBOX Test if required toolboxes are installed.
%
%   CHECKTOOLBOX tests whether toolboxes required by PRINCE are installed.
%
%   See also PRINCE


required_toobloxes = {'Statistics and machine learning toolbox'...
  'Curve fitting toolbox'...
  'Signal processing toolbox'...
  'Parallel computing toolbox'...
  'System identification toolbox'};

tmp = ver;
installed_toolboxes = {tmp(:).Name};

missing_toolboxes = required_toobloxes(~ismember(lower(required_toobloxes), ...
  lower(installed_toolboxes)));
if ~isempty(missing_toolboxes)
  ss = ['PrInCE needs the following Matlab toolboxes:\n\n    ' ...
    strjoin(missing_toolboxes, '\n    ') ...
    '\n\nWe can''t install them from here (sorry!). Please install them!'];
  error(sprintf(ss));
end


% reqtoolbox = [];
% 
% if ~license('test','statistics_toolbox')
%   reqtoolbox = [reqtoolbox '   - Statistics and machine learning toolbox\n'];
% end
% if ~license('test','curve_fitting_toolbox')
%   reqtoolbox = [reqtoolbox '   - Curve fitting toolbox\n'];
% end
% if ~license('test','signal_toolbox')
%   reqtoolbox = [reqtoolbox '   - Signal processing toolbox\n'];
% end
% if ~license('test','distrib_computing_toolbox')
%   reqtoolbox = [reqtoolbox '   - Parallel computing toolbox\n'];
% end
% if ~license('test','Identification_Toolbox')
%   reqtoolbox = [reqtoolbox '   - System identification toolbox\n'];
% end
% 
% % check toolboxes
% if ~isempty(reqtoolbox)
%   ss = ['PrInCE needs the following Matlab Toolboxes.\n' reqtoolbox '\n\nPlease install them!'];
%   error(sprintf(ss));
% end