clear all
% cd '/Users/nicholasoutram/Dropbox/Kernel estimation/Philippines/'

% I M P O R T   D A T A
[AD,AT] = xlsread('Data/data.xlsx', 'Adults'); % from excel file
[JD,JT] = xlsread('Data/data.xlsx','Juveniles');

DS = dir('GIS Polygons/Sampled Reefs/*shp'); % from sampled reefs
NSR = length(DS);
U = m_shaperead('GIS Polygons/Unsampled Reefs_WGS84_Pseudo'); % from unsampled reefs
NUR = length(U.ncst);

% F O R M A T   D A T A
    % REEFS DATA STRUCTURE: 
    % REEFNUM | LAT | LONG | ABUNDANCE | AREA | DENSITY
Reefs = AD(:,[1 3:5]);
Reefs(:,1) = int64(Reefs(:,1));
N = length(Reefs);

% S A M P L E D   R E E F S
AS = zeros(NSR,1);

for i = 1:NSR % calculate area of each sampled polygon
    SHPD = m_shaperead(['GIS Polygons/Sampled Reefs/' DS(i).name(1:end-4)]);
    Ix = SHPD.ncst{1}(:,1);
    Iy = SHPD.ncst{1}(:,2);
    
    [Reefs(i,3), Reefs(i,2), ~] = centroid(Ix,Iy);
    AS(i) = polyarea(Ix,Iy);
end

% A P O   I S L A N D
% this only works because Apo Island is the last polygon in the sampled set
apo = 0; f= 1; % calculate Apo Island area
for j = 1:length(Ix)
    if isnan(Ix(j))
        apo = apo + polyarea(Ix(f:j-1), Iy(f:j-1)); f = j+1;
    end
end

Ix(isnan(Ix)) = [];
Iy(isnan(Iy)) = [];
[Reefs(26,3), Reefs(26,2), ~] = centroid(Ix,Iy); % calculate Apo Island centroid

Reefs(:,5) = AS; % create area column
Reefs(26,5) = apo; % manually fill out apo island area

avDens = nanmean(Reefs(:,4) ./ (Reefs(:,5))); % calc average density
Reefs(23:25,4) = avDens .* Reefs(23:25,5); % fill out missing values
Reefs(:,6) = Reefs(:,4) ./ Reefs(:,5);

% U N S A M P L E D   R E E F S
UnsampledMatrix = zeros(NUR,4);
UnsampledMatrix(:,1) = 1+NSR:NUR+NSR;
UnsampledMatrix(:,6) = avDens; % estimate average density for unsampled reefs

for i = 1:NUR % iterate over nonsampled reefs
    Ix = U.ncst{i}(:,1);
    Iy = U.ncst{i}(:,2);
    
    % split by NaN parts
    area = 0; f= 1;
    for j = 1:length(Ix)
        if isnan(Ix(j))
            area = area + polyarea(Ix(f:j-1), Iy(f:j-1)); f = j+1;
        end
    end
    
    if max(isnan(Ix)) % set area if no NaNs originally
        UnsampledMatrix(i,5) = area;
    else % set area if NaNs
        UnsampledMatrix(i,5) = polyarea(Ix,Iy); % area

    end
            
    % remove NaN entries
    Ix(isnan(Ix)) = []; Iy(isnan(Iy)) = [];
    
    [UnsampledMatrix(i, 3), UnsampledMatrix(i,2), ~] = centroid(Ix,Iy);
end

UnsampledMatrix(:,4) = UnsampledMatrix(:,5) .* UnsampledMatrix(:,6); % estimate abundance for unsampled reefs

% G E N E R A T E   P A R E N T A G E   M A T R I X
S = [1:N; hist(JD(:,1),1:N)]; maxS = max(S(2,:)); % calculate parentage matrix dimensions
ParentageMatrix = repmat(-1, NSR, maxS); % If we're including the unsampled reefs, set these values automatically to -1
%ParentageMatrix = repmat(nan, NSR, maxS); % If we're NOT including the unsampled reefs, set these values automatically to NAN

sample = int64(JD(:,[1 5])); % format sample list
sn = sum(S(2,:)); % number of samples
currentReef = 1; j = 1;

for i = 1:sn % loop over all samples
   
    % if move onto a new reef, start again from 1st column
    if currentReef ~= sample(i,1)
        ParentageMatrix(currentReef, j:end) = NaN; % rest of row is nan (no samples)
        j = 1;
    end
    
    % if parent reef is known, include in parentage matrix
    if sample(i,2)
        ParentageMatrix(currentReef, j) = sample(i,2);
    end
        
    j = j + 1; currentReef = sample(i,1);
end 

ParentageMatrix([23 25],:) = NaN; % no samples from reef 23 or 25
ParentageMatrix(currentReef, j:end) = NaN; % rest of final reef sample
% ( nan == no data / -1 == unknown )

% M O D I F Y   F O R   A D U L T   S A M P L I N G   P R O P O R T I O N
ArtificialReefs_Adults = Reefs; % this could change for future data e.g. sample different reefs for adults / juveniles
SampledProp_Adults = ones(length(ArtificialReefs_Adults),1) .* .2; % this could be extracted from the data

if sum(SampledProp_Adults) > 0
    ArtificialReefs_Adults(:,5) = ArtificialReefs_Adults(:,5) .* (1 - SampledProp_Adults);
    Reefs(:,5) = Reefs(:,5) .* SampledProp_Adults;
    
% J O I N   M A T R I C E S

    AllReefs = vertcat(Reefs,UnsampledMatrix,ArtificialReefs_Adults);
else
    AllReefs = vertcat(Reefs, UnsampledMatrix);
end

% F O R M A T   A D D I T I O N A L   T H I N G S
reef_size = AllReefs(:,5) ./ 1000000; % convert area to km^2
Pop_dens = AllReefs(:,6);

% G E N E R A T E   D I S T A N C E   M A T R I X
D = squareform(pdist(AllReefs(:,[3 2]))) ./ 1000; % and convert to km

% S A V E   T O   F I L E
save philippinesData.mat ParentageMatrix D reef_size Pop_dens