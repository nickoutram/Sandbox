function [BestK] = fitDataP3(True_k, replicates)

if nargin == 0 True_k = 0.2; replicates = 5; end

BestK = zeros(replicates,1);

load philippinesData.mat % load parentageMatrix, D, Pop_dens & reef_size

% extract parameters from data
NumReefs = length(D);
SampleSize_Juveniles = length(ParentageMatrix);
SampledReefs_Juveniles = 1:size(ParentageMatrix,1);
SampledReefs_Adults = 1:size(ParentageMatrix,1);
NonSampledSet_Adults = setdiff(1:NumReefs,SampledReefs_Adults);

KVec = linspace(3,10,100);

for k = 1:length(KVec)
    LL_Crosssection(k) = sub_Likelihood(KVec(k),NumReefs,D,Pop_dens,reef_size,SampledReefs_Adults, ...
        SampledReefs_Juveniles,SampleSize_Juveniles,NonSampledSet_Adults,ParentageMatrix);
end

figure(1), clf, 
subplot(2,1,1), hold on, plot(KVec,LL_Crosssection)
[~,INDEX] = min(LL_Crosssection);
plot(KVec(INDEX),LL_Crosssection(INDEX),'ro')

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
BINS = [1:NumReefs];
k_f = exp(k_f);

% Construct a connectivity matrix based on the proposed kernel parameter
C = exp(-D.^2./k_f);

% Calculate the transport matrix that would be implied by the population sizes and connectivity matrix
L = repmat(Pop_dens.*reef_size,1,NumReefs).*C;

% WE GET ZEROS in L from C?? (couldn't be from Pop_dens or reef_size)

% Normalise the recruits (we're interested in proportions)
Prop = L./repmat(sum(L),NumReefs,1); % At each reef, this is the proportion of recruits from each reef

% If we didn't sample all reefs for adults, we have to worry about juveniles who can't be assigned to a source reef
L_mat = ones(length(SampledReefs_Juveniles),SampleSize_Juveniles);

% This loop is over the subset of destination reefs where juveniles were sampled
for i = 1:length(SampledReefs_Juveniles)
    % Nominate the sampled juvenile reefs one-by-one
    ThisJuvReef = SampledReefs_Juveniles(i);
    ThisRow = ParentageMatrix(ThisJuvReef,:);
    
    ThisRow(isnan(ThisRow)==1) = []; % remove NaN entries
    
    Unknown_likelihood = sum(Prop(NonSampledSet_Adults,ThisJuvReef)); % calc unknown likelihood - WHAT DOES IT RLY MEAN

    keyboard
    if isempty(ThisRow) % if no valid samples set multinomial to 1 to cancel out for this reef
        MN(i,1) = 1;
    
    else        
        % Calculate the number of samples for this reef by the number of non-NAN entries
        SampleSize_Juveniles_ThisReef = length(ThisRow);
        
        for j = 1:SampleSize_Juveniles_ThisReef % only loop through non NaN entries
            if ThisRow(j) == -1 % remove this clause if not including unknowns
                L_mat(i,j) = Unknown_likelihood;
                ThisRow(j) = 1; % change entry in ThisRow to NaN to remove for multinomial(?) LOOK AT DIS
            else
                L_mat(i,j) = Prop(ThisRow(j),ThisJuvReef);
            end
        end
%         ThisRow(isnan(ThisRow)==1) = [];
%         MN(i,1) = multinomial(length(ThisRow),histc(ThisRow,BINS));
        MN(i,1) = multinomial(SampleSize_Juveniles_ThisReef,histc(ThisRow,BINS));
    end
end
LL = -sum(log(MN.*prod(L_mat,2)));
