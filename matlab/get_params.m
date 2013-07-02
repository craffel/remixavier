function [P,X] = get_params(args, table, progname)
% [P,X] = get_params(args, table, progname)
%    Takes a cell array of name, default, description triples 
%    and parses args to set '-name value' args
%    or uses the default to return a P.name field in P.
%    for -help, or if args are unrecognized, generates a help
%    message. 
%    X returns any excess arguments (those not starting with -)
% 2013-07-01 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3; progname = 'get_params:'; end

warn = 1;

n = length(table);
if (mod(n, 3))
  error('Each option must be a string/value/description triple.');
end

P = struct();
% Set up all the defaults
for j = 1:3:n
  P = setfield(P,table{j},table{j+1});
end

X = {};

% Now process all arguments
nunused = 0;
i = 1;
while i <= length(args)
  found = 0;
  if args{i}(1) == '-'
    for j=1:3:n
      if strcmpi(args{i}(2:end), table{j})
        % dpwe: promote string args to numeric if default arg is
        % numeric
        if isnumeric(table{j+1}) && ischar(args{i+1})
          val = str2num(args{i + 1});
          % this occurs when specifying arguments from command line
        else
          % in all other cases, take what you get
          val = args{i + 1};
        end
        found = 1;
        P = setfield(P,table{j},val);
        break;
      end
    end
    if (~found)
      nunused = nunused + 1;
      if (warn)
        warning(sprintf('Option ''%s'' not used.', args{i}));
      else
        unused{2 * nunused - 1} = args{i};
        if length(args) > i
          % don't demand a value for the last, unused tag (e.g. -help)
          unused{2 * nunused} = args{i + 1};
        end 
      end
    end
    i = i+2;
  else
    % Arg that did not begin with '-'
    X{length(X)+1} = args{i};
    i = i+1;
  end
end

% Maybe print error message
if nunused > 0
  disp(progname);
  for j = 1:3:n
    if isnumeric(table{j+1})
      val = num2str(table{j+1});
    else
      val = table{j+1};
    end
    disp(sprintf('  -%s\t%s (%s)', table{j}, table{j+2}, val));
  end
  error('unrecognized arguments');
end
