function [BestK] = fitData(True_k, replicates)

if nargin == 0 True_k = 0.05; replicates = 5; end
BestK = zeros(replicates,1);
% true_k is only used for search bounds and graphing so may not be necessary.. check again later.

load philippinesData.mat % load parentageMatrix, D, Pop_dens & reef_size
NumReefs = length(D); SampleSize_Juveniles = length(ParentageMatrix);
SampledReefs_Juveniles = 1:size(ParentageMatrix,1); % do we need to check for NaN rows?
SampledReefs_Adults = 1:size(ParentageMatrix,1);
NonSampledSet_Adults = setdiff(1:NumReefs,SampledReefs_Adults);

KVec = linspace(3,10,100);

for k = 1:length(KVec)
    LL_CrossSection(k) = sub_Likelihood(KVec(k),NumReefs,D,Pop_dens,reef_size,SampledReefs_Adults, ...
        SampledReefs_Juveniles,SampleSize_Juveniles,NonSampledSet_Adults,ParentageMatrix);
end

figure(1), clf, subplot(2,1,1), hold on
plot(KVec,LL_CrossSection)
[~, INDEX] = min(LL_CrossSection);
plot(KVec(INDEX), LL_CrossSection(INDEX), 'ro')

subplot(2,1,2)
X = linspace(0,200,1000);
Y = exp(-X.^2./exp(KVec(INDEX)));
plot(X,Y,'r','linewidth',2)
return

for z = 1:replicates
    [Fit_k,~] = fminbnd(@sub_Likelihood,5*log(True_k),0.1*log(True_k),[], ...
        NumReefs,D,Pop_dens,reef_size,SampledReefs_Adults,SampledReefs_Juveniles,SampleSize_Juveniles,NonSampledSet_Adults,ParentageMatrix);
    
    % if fminbnd returns a search boundary change to NaN
    if Fit_k >= 1.01 * 0.1 * log(True_k) || Fit_k <= .99 * 5 * log(True_k)
       Fit_k = NaN;
    end
    
    BestK(z,1) = exp(Fit_k); % pass back actual value
end

function LL = sub_Likelihood(k_f,NumReefs,D,Pop_dens,reef_size,SampledReefs_Adults, ...
    SampledReefs_Juveniles,SampleSize_Juveniles,NonSampledSet_Adults,ParentageMatrix)

% if unknowns are included in sample include them in the multinomial?
if min(min(ParentageMatrix)) == -1
    BINS = [1:NumReefs]; % [-1 1:NumReefs]; 
else
    BINS = [1:NumReefs];
end

k_f = exp(k_f);
 
% Construct a connectivity matrix based on the proposed kernel parameter
C = exp(-D.^2./k_f); % Construct a connectivity matrix based on the proposed kernel parameter

% Calculate the transport matrix that would be implied by the population sizes and connectivity matrix
L = repmat(Pop_dens.*reef_size,1,NumReefs).*C;
 
% Calculate the transport matrix that would be implied by just the connectivity matrix
%  L = repmat(ones(size(reef_size)),1,NumReefs).*C; .. unsure what the point of this was
 
% Normalise the recruits (we're interested in proportions)
Prop = L./repmat(sum(L),NumReefs,1); % At each reef, this is the proportion of recruits from each reef
 
% If we didn't sample all reefs for adults, we have to worry about juveniles who can't be assigned to a source reef
L_mat = ones(length(SampledReefs_Juveniles),SampleSize_Juveniles);

% This loop is over the subset of destination reefs where juveniles were sampled
for i = 1:length(SampledReefs_Juveniles)
    ThisJuvReef = SampledReefs_Juveniles(i); % nominate the sampled juvenile reefs one-by-one
    thisReef_Samples = ParentageMatrix(ThisJuvReef,:);
    thisReef_Samples(isnan(thisReef_Samples)==1) = []; % remove NaN entries
    
    if isempty(thisReef_Samples) % if no valid samples, set multinomial to 1 to cancel out for reef
        MN(i,1) = 1;
        
    else
        % Calculate the number of samples for this reef by the number of non-NAN entries
        SampleSize_Juveniles_ThisReef = length(thisReef_Samples);
        
        for j = 1:SampleSize_Juveniles_ThisReef % only loop through non NaN entries
            if thisReef_Samples(j) == -1
                Unknown_likelihood = sum(Prop(NonSampledSet_Adults,ThisJuvReef)); % the likelihood that a juvenile came from a reef where adults weren't sampled
                L_mat(i,j) = Unknown_likelihood;
                % need to adjust for multinomial I think?
                %thisReef_Samples(j) = 1; % change entry in samples list to for multinomial? is this the correct number??
            else
                L_mat(i,j) = Prop(thisReef_Samples(j), ThisJuvReef);
            end
        end
        
        MN(i,1) = multinomial(SampleSize_Juveniles_ThisReef,histc(thisReef_Samples,BINS));
    end
end

LL = -sum(log(MN.*prod(L_mat,2)));