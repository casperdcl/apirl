%parpool('local', 6);
i=0;
for subj = [5, 20]
for dose = [0.1, 1, 10, -1, .7, 0.2, 0.4]
i = i + 1

if subj==5, tumours=-1;
else, tumours=-2;
end

% negative counts -> no noise
% negative tumours -> deterministic tumour generation
reconAPIRL = MLEM(1, 1.5, 0.75, sprintf('subject_%02d', subj), ...
                  430e6 * dose, tumours);

end, end
