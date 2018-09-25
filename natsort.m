function [X,ndx,dbg] = natsort(X,xpr,varargin) %#ok<*SPERR>
% Alphanumeric / Natural-Order sort the strings in a cell array of strings.
%
% (c) 2017 Stephen Cobeldick
%
% Alphanumeric sort of a cell array of strings: sorts by character order
% and also by the values of any numbers that are within the strings. The
% default is case-insensitive ascending with integer number substrings:
% optional inputs control the sort direction, case sensitivity, and the
% number substring matching (see the section "Number Substrings" below).
%
%%% Syntax:
%  Y = natsort(X)
%  Y = natsort(X,xpr)
%  Y = natsort(X,xpr,<options>)
% [Y,ndx] = natsort(X,...);
%
% To sort filenames or filepaths use NATSORTFILES (File Exchange 47434).
% To sort the rows of a cell array of strings use NATSORTROWS (File Exchange 47433).
%
% See also NATSORTFILES NATSORTROWS SORTROWS SORT CELLSTR REGEXP SSCANF NUM2ORDINAL NUM2WORDS NUM2BIP NUM2SIP INTMAX
%
%% Number Substrings %%
%
% By default consecutive digit characters are interpreted as an integer.
% The optional regular expression pattern <xpr> permits the numbers to also
% include a +/- sign, a decimal point, exponent E-notation or any literal
% characters, quantifiers or look-around requirements. For more information:
% http://www.mathworks.com/help/matlab/matlab_prog/regular-expressions.html
%
% The substrings are then parsed by SSCANF into numeric variables, using
% either the *default format '%f', or the user-supplied format specifier.
%
% This table shows some example regular expression patterns for some common
% notations and ways of writing numbers (see section "Examples" for more):
%
% <xpr> Regular | Number Substring | Number Substring              | SSCANF
% Expression:   | Match Examples:  | Match Description:            | Format Specifier:
% ==============|==================|===============================|==================
% *         \d+ | 0, 1, 234, 56789 | unsigned integer              | %f  %u  %lu  %i
% --------------|------------------|-------------------------------|------------------
%     (-|+)?\d+ | -1, 23, +45, 678 | integer with optional +/- sign| %f  %d  %ld  %i
% --------------|------------------|-------------------------------|------------------
%   \d+(\.\d+)? | 012, 3.45, 678.9 | integer or decimal            | %f
% --------------|------------------|-------------------------------|------------------
%   \d+|Inf|NaN | 123, 4, Inf, NaN | integer, infinite or NaN value| %f
% --------------|------------------|-------------------------------|------------------
%  \d+\.\d+e\d+ | 0.123e4, 5.67e08 | exponential notation          | %f
% --------------|------------------|-------------------------------|------------------
%  0[0-7]+      | 012, 03456, 0700 | octal prefix & notation       | %o  %i
% --------------|------------------|-------------------------------|------------------
%  0X[0-9A-F]+  | 0X0, 0XFF, 0X7C4 | hexadecimal prefix & notation | %x  %i
% --------------|------------------|-------------------------------|------------------
%  0B[01]+      | 0B101, 0B0010111 | binary prefix & notation      | %b
% --------------|------------------|-------------------------------|------------------
%
% The SSCANF format specifier (including %b) can include literal characters
% and skipped fields. The octal, hexadecimal and binary prefixes are optional.
% For more information: http://www.mathworks.com/help/matlab/ref/sscanf.html
%
%% Relative Sort Order %%
%
% The sort order of the number substrings relative to the characters
% can be controlled by providing one of the following string options:
%
% Option Token:| Relative Sort Order:                 | Example:
% =============|======================================|====================
% 'beforechar' | numbers < char(0:end)                | '1' < '.' < 'A'
% -------------|--------------------------------------|--------------------
% 'afterchar'  | char(0:end) < numbers                | '.' < 'A' < '1'
% -------------|--------------------------------------|--------------------
% 'asdigit'   *| char(0:47) < numbers < char(48:end)  | '.' < '1' < 'A'
% -------------|--------------------------------------|--------------------
%
% Note that the digit characters have character values 48 to 57, inclusive.
%
%% Examples %%
%
%%% Integer number substrings:
% A = {'a2', 'a10', 'a1'};
% sort(A)
%  ans =  'a1'  'a10'  'a2'
% natsort(A)
%  ans =  'a1'  'a2'  'a10'
%
%%% Multiple number substrings (e.g. release version numbers):
% B = {'v10.6', 'v9.10', 'v9.5', 'v10.10', 'v9.10.20', 'v9.10.8'};
% sort(B)
%  ans =  'v10.10'  'v10.6'  'v9.10'  'v9.10.20'  'v9.10.8'  'v9.5'
% natsort(B)
%  ans =  'v9.5'  'v9.10'  'v9.10.8'  'v9.10.20'  'v10.6'  'v10.10'
%
%%% Integer, decimal or Inf number substrings, possibly with +/- signs:
% C = {'test+Inf', 'test11.5', 'test-1.4', 'test', 'test-Inf', 'test+0.3'};
% sort(C)
%  ans =  'test'  'test+0.3'  'test+Inf'  'test-1.4'  'test-Inf'  'test11.5'
% natsort(C, '(-|+)?(Inf|\d+(\.\d+)?)')
%  ans =  'test'  'test-Inf'  'test-1.4'  'test+0.3'  'test11.5'  'test+Inf'
%
%%% Integer or decimal number substrings, possibly with an exponent:
% D = {'0.56e007', '', '4.3E-2', '10000', '9.8'};
% sort(D)
%  ans =  ''  '0.56e007'  '10000'  '4.3E-2'  '9.8'
% natsort(D, '\d+(\.\d+)?(E(+|-)?\d+)?')
%  ans =  ''  '4.3E-2'  '9.8'  '10000'  '0.56e007'
%
%%% Hexadecimal number substrings (possibly with '0X' prefix):
% E = {'a0X7C4z', 'a0X5z', 'a0X18z', 'aFz'};
% sort(E)
%  ans =  'a0X18z'  'a0X5z'  'a0X7C4z'  'aFz'
% natsort(E, '(?<=a)(0X)?[0-9A-F]+', '%x')
%  ans =  'a0X5z'  'aFz'  'a0X18z'  'a0X7C4z'
%
%%% Binary number substrings (possibly with '0B' prefix):
% F = {'a11111000100z', 'a0B101z', 'a0B000000000011000z', 'a1111z'};
% sort(F)
%  ans =  'a0B000000000011000z'  'a0B101z'  'a11111000100z'  'a1111z'
% natsort(F, '(0B)?[01]+', '%b')
%  ans =  'a0B101z'  'a1111z'  'a0B000000000011000z'  'a11111000100z'
%
%%% uint64 number substrings (with full precision!):
% natsort({'a18446744073709551615z', 'a18446744073709551614z'}, [], '%lu')
%  ans =   'a18446744073709551614z'  'a18446744073709551615z'
%
%%% Case sensitivity:
% G = {'a2', 'A20', 'A1', 'a10', 'A2', 'a1'};
% natsort(G, [], 'ignorecase') % default
%  ans =   'A1'  'a1'  'a2'  'A2'  'a10'  'A20'
% natsort(G, [], 'matchcase')
%  ans =   'A1'  'A2'  'A20'  'a1'  'a2'  'a10'
%
%%% Sort direction:
% H = {'2', 'a', '3', 'B', '1'};
% natsort(H, [], 'ascend') % default
%  ans =   '1'  '2'  '3'  'a'  'B'
% natsort(H, [], 'descend')
%  ans =   'B'  'a'  '3'  '2'  '1'
%
%%% Relative sort-order of number substrings compared to characters:
% X = num2cell(char(32+randperm(63)));
% cell2mat(natsort(X, [], 'asdigit')) % default
%  ans = '!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_'
% cell2mat(natsort(X, [], 'beforechar'))
%  ans = '0123456789!"#$%&'()*+,-./:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_'
% cell2mat(natsort(X, [], 'afterchar'))
%  ans = '!"#$%&'()*+,-./:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_0123456789'
%
%% Input and Output Arguments %%
%
%%% Inputs (*=default):
%   X   = Cell of Strings, with strings to be sorted into natural-order.
%   xpr = String Token, regular expression to detect number substrings, *'\d+'.
% <options> string tokens can be entered in any order, as many as required:
%   - Case sensitive/insensitive matching: 'matchcase'/'ignorecase'*.
%   - Sort direction: 'descend'/'ascend'*.
%   - Relative sort of numbers: 'beforechar'/'afterchar'/'asdigit'*.
%   - The SSCANF number conversion format, e.g.: *'%f', '%x', '%i', etc.
%
%%% Outputs:
%   Y   = Cell of Strings, <X> with all strings sorted into natural-order.
%   ndx = Numeric Array, such that Y = X(ndx). The same size as <X>.
%   dbg = Cell Array of all parsed characters and number values. Each row is
%         one string, linear-indexed from <X>. Helps to debug string parsing.
%
% [X,ndx,dbg] = natsort(X,*xpr,<options>)

%% Input Wrangling %%
%
assert(iscell(X),'First input <X> must be a cell array.')
tmp = cellfun('isclass',X,'char') & 2>cellfun('size',X,1) & 3>cellfun('ndims',X);
assert(all(tmp(:)),'First input <X> must be a cell array of strings (1xN character).')
%
% Regular expression:
if nargin<2 || isnumeric(xpr)&&isempty(xpr)
	xpr = '\d+';
else
	assert(ischar(xpr)&&isrow(xpr),'Second input <xpr> must be a regular expression.')
end
%
% Optional arguments:
tmp = cellfun('isclass',varargin,'char') & 1==cellfun('size',varargin,1) & 2==cellfun('ndims',varargin);
assert(all(tmp(:)),'All optional arguments must be strings (1xN character).')
% Character case matching:
MatL = strcmpi(varargin,'matchcase');
CasL = strcmpi(varargin,'ignorecase')|MatL;
% Sort direction:
DesL = strcmpi(varargin,'descend');
DirL = strcmpi(varargin,'ascend')|DesL;
% Relative sort-order of numbers compared to characters:
BefL = strcmpi(varargin,'beforechar');
AftL = strcmpi(varargin,'afterchar');
RsoL = strcmpi(varargin,'asdigit')|BefL|AftL;
% SSCANF conversion format:
FmtL = ~(CasL|DirL|RsoL);
%
if nnz(DirL)>1
	error('Sort direction is overspecified:%s\b.',sprintf(' ''%s'',',varargin{DirL}))
end
%
if nnz(RsoL)>1
	error('Relative sort-order is overspecified:%s\b.',sprintf(' ''%s'',',varargin{RsoL}))
end
%
FmtN = nnz(FmtL);
if FmtN>1
	error('Overspecified optional arguments:%s\b.',sprintf(' ''%s'',',varargin{FmtL}))
end
%
%% Split Strings %%
%
% Split strings into number and remaining substrings:
[MtS,MtE,MtC,SpC] = regexpi(X(:),xpr,'start','end','match','split',varargin{CasL});
%
% Determine lengths:
MtcD = cellfun(@minus,MtE,MtS,'UniformOutput',false);
LenZ = cellfun('length',X(:))-cellfun(@sum,MtcD);
LenY = max(LenZ);
LenX = numel(MtC);
%
dbg = cell(LenX,LenY);
NuI = false(LenX,LenY);
ChI = false(LenX,LenY);
ChA = char(double(ChI));
%
ndx = 1:LenX;
for k = ndx(LenZ>0)
	% Determine indices of numbers and characters:
	ChI(k,1:LenZ(k)) = true;
	if ~isempty(MtS{k})
		tmp = MtE{k} - cumsum(MtcD{k});
		dbg(k,tmp) = MtC{k};
		NuI(k,tmp) = true;
		ChI(k,tmp) = false;
	end
	% Transfer characters into char array:
	if any(ChI(k,:))
		tmp = SpC{k};
		ChA(k,ChI(k,:)) = [tmp{:}];
	end
end
%
%% Convert Number Substrings %%
%
if FmtN % One format specifier
	fmt = varargin{FmtL};
	err = ['Format specifier results in an empty output from sscanf: ''',fmt,''''];
	P = '(?<!%)(%%)*%'; % match an odd number of % characters.
	[T,S] = regexp(fmt,[P,'(\d*)(b|d|i|u|o|x|f|e|g|l(d|i|u|o|x))'],'tokens','split');
	assert(isscalar(T),'Unsupported optional argument: ''%s''',fmt)
	assert(isempty(T{1}{2}),'Format specifier cannot include field-width: ''%s''',fmt)
	switch T{1}{3}(1)
		case 'b' % binary
			fmt = regexprep(fmt,[P,'(\*?)b'],'$1%$2[01]');
			val = dbg(NuI);
			if numel(S{1})<2 || ~strcmpi('0B',S{1}(end-1:end))
				% Remove '0B' if not specified in the format string:
				val = regexprep(val,'(0B)?([01]+)','$2','ignorecase');
			end
			val = cellfun(@(s)sscanf(s,fmt),val, 'UniformOutput',false);
			assert(~any(cellfun('isempty',val)),err)
			NuA(NuI) = cellfun(@(s)sum(pow2(s-48,numel(s)-1:-1:0)),val);
		case 'l' % 64-bit
			NuA(NuI) = cellfun(@(s)sscanf(s,fmt),dbg(NuI)); %slow!
		otherwise % double
			NuA(NuI) = sscanf(sprintf('%s\v',dbg{NuI}),[fmt,'\v']); % fast!
	end
else % No format specifier -> double
	NuA(NuI) = sscanf(sprintf('%s\v',dbg{NuI}),'%f\v');
end
% Note: NuA's class is determined by SSCANF.
NuA(~NuI) = 0;
NuA = reshape(NuA,LenX,LenY);
%
%% Debugging Array %%
%
if nargout>2
	for k = reshape(find(NuI),1,[])
		dbg{k} = NuA(k);
	end
	for k = reshape(find(ChI),1,[])
		dbg{k} = ChA(k);
	end
end
%
%% Sort %%
%
if ~any(MatL) % ignorecase
	ChA = upper(ChA);
end
%
ide = ndx.';
% From the last column to the first...
for n = LenY:-1:1
	% ...sort the characters and number values:
	[C,idc] = sort(ChA(ndx,n),1,varargin{DirL});
	[~,idn] = sort(NuA(ndx,n),1,varargin{DirL});
	% ...keep only relevant indices:
	jdc = ChI(ndx(idc),n); % character
	jdn = NuI(ndx(idn),n); % number
	jde = ~ChI(ndx,n)&~NuI(ndx,n); % empty
	% ...define the sort-order of numbers and characters:
	jdo = any(AftL)|(~any(BefL)&C<48);
	% ...then combine these indices in the requested direction:
	if any(DesL) % descending
		idx = [idc(jdc&~jdo);idn(jdn);idc(jdc&jdo);ide(jde)];
	else % ascending
		idx = [ide(jde);idc(jdc&jdo);idn(jdn);idc(jdc&~jdo)];
	end
	ndx = ndx(idx);
end
%
ndx  = reshape(ndx,size(X));
X = X(ndx);
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%natsort