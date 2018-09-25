function [X,ndx] = natsortfiles(X,varargin)
% Alphanumeric / Natural-Order sort of a cell array of filenames/filepaths.
%
% (c) 2017 Stephen Cobeldick
%
% Alphanumeric sort of a cell array of filenames or filepaths: sorts by
% character order and also by the values of any numbers that are within
% the strings. Filenames, file-extensions, and directories (if supplied)
% are split apart and are sorted separately: this ensures that shorter
% filenames sort before longer ones (i.e. just like a dictionary sort).
%
%%% Syntax:
%  Y = natsortfiles(X)
%  Y = natsortfiles(X,xpr)
%  Y = natsortfiles(X,xpr,<options>)
% [Y,ndx] = natsortfiles(X,...)
%
% To sort all of the strings in a cell array use NATSORT (File Exchange 34464).
% To sort the rows of a cell array of strings use NATSORTROWS (File Exchange 47433).
%
% See also NATSORT NATSORTROWS SORT SORTROWS NUM2SIP NUM2BIP CELLSTR REGEXP DIR FILEPARTS FULLFILE FILESEP
%
%% File Dependency %%
%
% NATSORTFILES requires the function NATSORT (File Exchange 34464). The inputs
% <xpr> and <options> are passed directly to NATSORT: see NATSORT for case
% sensitivity, sort direction, numeric substring matching, and other options.
%
%% Explanation %%
%
% Using SORT on filenames will sort any of char(0:45), including the printing
% characters ' !"#$%&''()*+,-', before the file extension separator character '.'.
% Therefore this function splits the name and extension and sorts them separately.
%
% Similarly the file separator character within filepaths can cause longer
% directory names to sort before shorter ones, as char(0:46)<'/' and char(0:91)<'\'.
% NATSORTFILES splits filepaths at each file separator character and sorts
% every level of the directory hierarchy separately, ensuring that shorter
% directory names sort before longer, regardless of the characters in the names.
%
%% Examples %%
%
%%% DIR and Cell Array:
% D = 'C:\Test';
% S = dir(fullfile(D,'*.txt'));
% N = natsortfiles({S.name});
% for k = 1:numel(N)
%     fullfile(D,N{k})
% end
%
% A = {'a2.txt', 'a10.txt', 'a1.txt'};
% sort(A)         % Wrong numeric order:
%  ans = 'a1.txt'  'a10.txt'  'a2.txt'
% natsortfiles(A)
%  ans = 'a1.txt'  'a2.txt'  'a10.txt'
%
% B = {'test_new.m'; 'test-old.m'; 'test.m'};
% sort(B)         % Note '-' sorts before '.':
%  ans =
%    'test-old.m'
%    'test.m'
%    'test_new.m'
% natsortfiles(B) % Shorter names before longer (dictionary sort):
%  ans =
%    'test.m'
%    'test-old.m'
%    'test_new.m'
%
% C = {'test2.m'; 'test10-old.m'; 'test.m'; 'test10.m'; 'test1.m'};
% sort(C)         % Wrong numeric order:
%  ans =
%    'test.m'
%    'test1.m'
%    'test10-old.m'
%    'test10.m'
%    'test2.m'
% natsortfiles(C) % Correct numeric order, shorter names before longer:
%  ans =
%    'test.m'
%    'test1.m'
%    'test2.m'
%    'test10.m'
%    'test10-old.m'
%
%%% Directory Names:
% D = {'A2-old\test.m';'A10\test.m';'A2\test.m';'A1archive.zip';'A1\test.m'};
% sort(D)         % Wrong numeric order, and '-' sorts before '\':
%  ans =
%    'A10\test.m'
%    'A1\test.m'
%    'A1archive.zip'
%    'A2-old\test.m'
%    'A2\test.m'
% natsortfiles(D) % Shorter names before longer (dictionary sort):
%  ans =
%    'A1archive.zip'
%    'A1\test.m'
%    'A2\test.m'
%    'A2-old\test.m'
%    'A10\test.m'
%
%% Input and Output Arguments %%
%
% Please see NATSORT for a full description of <xpr> and the <options>.
%
%%% Inputs (*=default):
%  X   = Cell of Strings, with filenames or filepaths to be sorted.
%  xpr = String Token, regular expression to detect numeric substrings, '\d+'*.
%  <options> can be supplied in any order and are passed directly to NATSORT.
%
%%% Outputs:
%  Y   = Cell of Strings, <X> with the filenames sorted according to <options>.
%  ndx = Numeric Matrix, same size as <X>. Indices such that Y = X(ndx).
%
% [Y,ndx] = natsortrows(X,*xpr,<options>)

%% Input Wrangling %%
%
assert(iscell(X),'First input <X> must be a cell array.')
tmp = cellfun('isclass',X,'char') & cellfun('size',X,1)<2 & cellfun('ndims',X)<3;
assert(all(tmp(:)),'First input <X> must be a cell array of strings (1xN character).')
%
%% Split and Sort File Names/Paths %%
%
% Split full filepaths into file [path,name,extension]:
[pth,nam,ext] = cellfun(@fileparts,X(:),'UniformOutput',false);
% Split path into {dir,subdir,subsubdir,...}:
%pth = regexp(pth,filesep,'split'); % OS dependent
pth = regexp(pth,'[/\\]','split'); % either / or \
len = cellfun('length',pth);
vec(1:numel(len)) = {''};
%
% Natural-order sort of the file extension, file name, and directory names:
[~,ndx] = natsort(ext,varargin{:});
[~,ind] = natsort(nam(ndx),varargin{:});
ndx = ndx(ind);
for k = max(len):-1:1
	idx = len>=k;
	vec(~idx) = {''};
	vec(idx) = cellfun(@(c)c(k),pth(idx));
	[~,ind] = natsort(vec(ndx),varargin{:});
	ndx = ndx(ind);
end
%
% Return the sorted array and indices:
ndx = reshape(ndx,size(X));
X = X(ndx);
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%natsortfiles