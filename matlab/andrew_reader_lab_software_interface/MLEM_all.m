%parpool('local', 6);
i=0;
for dose = [0.1, 0.2, 0.4, 1, 10, -1]
for subj = [4, 54]
i=i+1, if i<=7, continue, end

if subj==4, tumours=-1;
else, tumours=-2;
end

reconAPIRL = MLEM(1, 1.5, 0.75, sprintf('subject_%02d', subj), ...
                  500e6 * dose, tumours);

end, end

