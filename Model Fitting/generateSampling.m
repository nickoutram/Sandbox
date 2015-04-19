function [samplingMatrix] = generateSampling(p, R)

% Generate samples from theorised distributions with kernel based larval
% dispersal
if nargin==0, p = 8.4444; R = 1; end

% generate transport matrix from parameters, and reef stats.
% do we want to include artificial reefs??? # YES?

load philippinesData.mat % load parentageMatrix, D, Pop_dens & reef_size
NumReefs = length(D);
SampledReefs_Juveniles = 1:size(ParentageMatrix,1);
BINS = [-1 SampledReefs_Juveniles];

k_f = exp(p); % parameterise kernel parameter
C = exp(-D.^2 ./ k_f); % connectivity matrix
L = repmat(Pop_dens .* reef_size, 1, NumReefs) .* C; % transport matrix
Prop = L ./ repmat(sum(L), NumReefs, 1); % normalised proportion matrix for probabilities

% determine the number of samples taken on each reef
ReefSamples = zeros(length(SampledReefs_Juveniles),1);
for thisReef = 1:length(SampledReefs_Juveniles)
    thisReef_Samples = ParentageMatrix(thisReef,:);
    thisReef_Samples(isnan(thisReef_Samples) == 1) = [];
    ReefSamples(thisReef) = length(thisReef_Samples);
end

% random sample with replacement in each reef
SampleMatrix = NaN(length(SampledReefs_Juveniles), size(ParentageMatrix,2), R); % empty sample matrix

for r = 1:R % generate R sets of samples over SampledReefs_Juveniles reefs
    for s = 1:length(SampledReefs_Juveniles)
        SampleMatrix(s, 1:ReefSamples(s), r) = randsample(NumReefs, ReefSamples(s), true, Prop(:,s));
    end
end

SampleMatrix(SampleMatrix > length(SampledReefs_Juveniles)) = -1; % replace unknown reefs with -1
sumSampleMatrix = permute(histc(permute(SampleMatrix, [2 1 3]), BINS), [2 1 3]); % sum number from each reef (sampled at each reef)
samplingMatrix = sumSampleMatrix;
end