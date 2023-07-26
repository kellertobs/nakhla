function progressbar(input, varargin)
%PROGRESSBAR Display a progress bar inside the command window
%   64%  [##############################..................]
% 
%   Downloaded from https://github.com/YarmoM/progress.git
%   19 April 2022
%
%   progressbar('_start') initializes a new progress bar. Must always be
%   called first.
%   
%   progressbar(i) updates the progress bar. i is a percentage.
%   progressbar(i, m) is similar, but a percentage is automagically
%   calculated, where i is the current step and m the maximum number of
%   steps.
%   
%   progressbar('_end') ends the progress bar.
%   
%   progressbar('_erase') ends the progress bar and removes if from the
%   command window, restoring it to the way it was before the progress
%   bar was initiated.
%   
%   progressbar(message) ends the progress bar and adds message after the bar.
%
%   progressbar(..., opts) uses options from the opts struct. See below.
%
%   Option settings:
%       opts.percentageLength: sets the number of characters reserved for
%           the percentage display (default: 5)
%       opts.barLength: sets the number of characters reserved for the
%           progress bar (default: 48)
%       opts.charEmpty: sets the "empty" character (default: '.')
%       opts.charFilled: sets the "filled" character (default: '#')

    % Handle input
    p = inputParser;
    p.addRequired('input');
    p.addOptional('max', []);
    p.addOptional('opts', []);
    p.parse(input, varargin{:});
    p = p.Results;

    % Handle options
    opts = struct;
    opts.percentageLength = defaultOption(p.opts, 'percentageLength', 5);
    opts.barLength = defaultOption(p.opts, 'barLength', 20);
    opts.charEmpty = defaultOption(p.opts, 'charEmpty', '.');
    opts.charFilled = defaultOption(p.opts, 'charFilled', '#');

    % Act depending on the input
    if ischar(p.input) % If input is a string
        switch p.input
            case '_start' % Start a new progress bar
                strOut = generateString(0, opts);
                fprintf(strOut);
            case '_end' % End a progress bar
            case '_erase' % End and erase a progress bar
                strOut = generateString(0, opts);
                strCR = repmat('\b',1,length(strOut)-2);
                fprintf(strCR);
            otherwise % End a progress bar with a message after the bar
                strOut = sprintf(' %s\n', p.input);
                strCR = repmat('\b',1,1);
                fprintf([strCR strOut]);
        end
    elseif isnumeric(p.input) % If input is a number
        % If a max value is also provided, compute the percentage
        % If no, p.input is assumed to be a percentage
        if ~isempty(p.max) && isnumeric(p.max)
            p.input = round(100*input/p.max);
        else
            if p.input < 0
                p.input = 0;
            elseif p.input > 100
                p.input = 100;
            end
        end
        
        % Generate the output string
        strOut = generateString(p.input, opts);
        strCR = repmat('\b',1,length(strOut)-2);
        fprintf([strCR strOut]);
    else
        % Unsupported argument type
        error('Unsupported argument type');
    end
end

function str = generateString(percentage, s)
% Generate the progress bar string

    percentage = floor(percentage);
    percentageOut = [num2str(percentage) '%%'];
    percentageOut = [percentageOut repmat(' ',1,s.percentageLength-length(percentageOut)+1)];
    nDots = floor(percentage/100*s.barLength);
    dotOut = ['[' repmat('#',1,nDots) repmat('.',1,s.barLength-nDots) ']'];
    str = [percentageOut dotOut '\n'];
end

function output = defaultOption(opts, field, default)
% Check whether a setting was found in the options struct

    if isfield(opts, field)
        output = opts.(field);
    else
        output = default;
    end
end