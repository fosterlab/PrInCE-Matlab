function reqtoolbox = checktoolbox()

% CHECKTOOLBOX Test if required toolboxes are installed.
%
%   CHECKTOOLBOX tests whether toolboxes required by PRINCE are installed.
%
%   See also PRINCE


reqtoolbox = [];

if ~license('test','statistics_toolbox')
  reqtoolbox = [reqtoolbox '   - Statistics and machine learning toolbox\n'];
end
if ~license('test','curve_fitting_toolbox')
  reqtoolbox = [reqtoolbox '   - Curve fitting toolbox\n'];
end
if ~license('test','signal_toolbox')
  reqtoolbox = [reqtoolbox '   - Signal processing toolbox\n'];
end
if ~license('test','distrib_com2puting_toolbox')
  reqtoolbox = [reqtoolbox '   - Parallel computing toolbox\n'];
end
if ~license('test','Identification_Toolbox')
  reqtoolbox = [reqtoolbox '   - System identification toolbox\n'];
end

% check toolboxes
if ~isempty(reqtoolbox)
  ss = ['PrInCE needs the following Matlab Toolboxes.\n' reqtoolbox '\n\nPlease install them!'];
  error(sprintf(ss));
end