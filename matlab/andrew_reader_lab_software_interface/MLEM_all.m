%%
%{
clear all, close all, delete(gcp('nocreate'))
parpool('local', 6);
%}
%%

i=0;
%for tumours = [0 1]
for tumours = [0]
for subj = [1, 2]
%for subj = [5, 20]
%for dose = [0.1, 1, 10, -1, .7, 0.2, 0.4]
for dose = [0.1, 0.7]
i = i + 1
%if i < 3, continue, end

%if subj==5, tumours=tumours*-1;
%else, tumours=tumours*-2;
%end

% negative counts -> no noise
% negative tumours -> deterministic tumour generation

%% brainweb reconstruction simulation
%reconAPIRL = MLEM(1, 1, 0.75, sprintf('subject_%02d', subj), ...
%  430e6 * dose, tumours, ...
%  '/scratch/cc16/apirl/bw/0', [4 6 7] + 1);
%% real patient re-reconstruction
reconAPIRL = MLEM(1, 1, 0.75, sprintf('AD_%d', subj), ...
  430e6 * dose, -subj * tumours, ...
  '/scratch/cc16/apirl/0', [4 6 7] + 1);

end, end
end
