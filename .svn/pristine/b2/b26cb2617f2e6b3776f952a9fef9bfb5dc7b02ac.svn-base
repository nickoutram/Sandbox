function [BestK] = fitDataP(True_k, p, replicates)

%unpack parameters
p = num2cell(p);
[SampleSize_Juveniles, PropSampledJuvenileReefs, NumAdultReefs, SampleSize_Adults] = p{:};
BestK = zeros(replicates,1);

load philippinesData.mat

for z = 1:replicates
    % Get parameters from the output of the populkation model
    NumReefs = length(D);
    
    NumberSampledJuvenileReefs = size(ParentageMatrix,1);
    SampledReefs_Juveniles = % count non NaN rows in ParentageMatrix
    
    SampledReefs_Adults = size(ParentageMatrix,1);
    
    NonSampledSet_Adults = setdiff(1:NumReefs,SampledReefs_Adults);

    [Fit_k,~] = fminbnd(@sub_Likelihood,5*log(True_k),0.1*log(True_k),[], ...
        NumReefs,D,Pop_dens,reef_size,SampledReefs_Adults,SampledReefs_Juveniles,SampleSize_Juveniles,NonSampledSet_Adults,ParentageMatrix);
    
    % if fminbnd returns a search boundary change to NaN
    if Fit_k >= 1.01 * 0.1 * log(True_k) || Fit_k <= .99 * 5 * log(True_k)
       Fit_k = NaN;
    end
    
    % pass back actual value
    BestK(z,1) = exp(Fit_k); 
    
%     s = [];
%     for k_f = linspace(10*log(True_k),0.1.*log(True_k),100)
%         s = [s sub_Likelihood(k_f,NumReefs,D,Pop_dens,reef_size,SampledReefs_Adults,SampledReefs_Juveniles,SampleSize_Juveniles,NonSampledSet_Adults,ParentageMatrix)];
%     end, plot(linspace(10*log(True_k),0.1.*log(True_k),100),s,'r','linewidth',2)
%     s = [];
%     for k_f = linspace(5*log(True_k),0.1.*log(True_k),50)
%         s = [s sub_Likelihood(k_f,NumReefs,D,Pop_dens,reef_size,SampledReefs_Adults,SampledReefs_Juveniles,SampleSize_Juveniles,NonSampledSet_Adults,ParentageMatrix)];
%     end, plot(linspace(5*log(True_k),0.1.*log(True_k),50),s,'b--','linewidth',2)
%     pause, cla
end


function LL = sub_Likelihood(k_f,NumReefs,D,Pop_dens,reef_size,SampledReefs_Adults, ...
    SampledReefs_Juveniles,SampleSize_Juveniles,NonSampledSet_Adults,ParentageMatrix)
BINS = 1:NumReefs;
k_f = exp(k_f);
 
% Construct a connectivity matrix based on the proposed kernel parameter
C = exp(-D.^2./k_f);

% Calculate the transport matrix that would be implied by the population sizes and connectivity matrix
L = repmat(Pop_dens.*reef_size,1,NumReefs).*C;
 
% Calculate the transport matrix that would be implied by just the connectivity matrix
%     L = repmat(ones(size(reef_size)),1,NumReefs).*C;
 
% Normalise the recruits (we're interested in proportions)
Prop = L./repmat(sum(L),NumReefs,1); % At each reef, this is the proportion of recruits from each reef
 
% If we didn't sample all reefs for adults, we have to worry about juveniles
% who can't be assigned to a source reef
L_mat = zeros(length(SampledReefs_Juveniles),SampleSize_Juveniles);

% This loop is over the subset of destination reefs where juveniles were sampled
for i = 1:length(SampledReefs_Juveniles)
    % Nominate the sampled juvenile reefs one-by-one
    ThisJuvReef = SampledReefs_Juveniles(i);
    
    % This is the likelihood that a juvenile came from a reef where adults weren't sampled
    Unknown_likelihood = sum(Prop(NonSampledSet_Adults,ThisJuvReef));

    % Calculate the number of samples for this reef by the number of non-NAN entries
    SampleSize_Juveniles_ThisReef = sum(isnan(ParentageMatrix(ThisJuvReef,:)==0));
    
    % This loop is over the samples
    for j = 1:SampleSize_Juveniles_ThisReef
        if isnan(ParentageMatrix(ThisJuvReef,j)) == 0
            if ismember(ParentageMatrix(ThisJuvReef,j),SampledReefs_Adults) == 1
                L_mat(i,j) = Prop(ParentageMatrix(ThisJuvReef,j),ThisJuvReef);
            else
                L_mat(i,j) = Unknown_likelihood;
            end
        else L_mat(i,j) = 1;
        end
    end
    if sum(ParentageMatrix(ThisJuvReef,:)) == 0
        MN(i,1) = 1;
        L_mat(i,:) = 1;
    else
        MN(i,1) = multinomial(SampleSize_Juveniles_ThisReef,histc(ParentageMatrix(ThisJuvReef,:),BINS));
    end
end
LL = -sum(log(MN.*prod(L_mat,2)));







