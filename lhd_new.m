% VERSION180804
%
function x_start = lhd_new(size_pop,vlb,vub)

% --------------------------------------------------------------------------
% parameter explanation and settings
% --------------------------------------------------------------------------
% input parameters:
% size_pop    population size (scalar)
% vlb         lower search bounds (array)
% vub         upper search bounds (array)

% output parameter
% x_start     a collection of vectors in the parameter space

% LHD stands for Latin Hypercube Design
% made by Camiel on 25 sept 2002
% LHD creates a random series of starting vectors,
% that are well spread over the parameter space

% new version made on 25-2-2004
% remark that the created population does NOT satisfy extra requirements!

stap = (vub - vlb)/size_pop;
x = ones(size_pop,1)*vlb;

dx = (1:size_pop)'*stap;
disturb = rand(size_pop,length(vlb)).*(ones(size_pop,1)*stap);
dx = dx - disturb;

x = x + dx;
for param = 1:length(vlb)
    ind = randperm(size_pop)';
    x(:,param) = x(ind,param);
end
x_start = x;
