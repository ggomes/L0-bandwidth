% script to load and plot data from Bullock/Day study.

clear
close all
add_dependencies()

rootfolder = fileparts( mfilename('fullpath'));
datafolder = fullfile(rootfolder,'noblesville_data');
matlabdatafile.G = fullfile(datafolder,'G');
matlabdatafile.veharrivals_agg = fullfile(datafolder,'veharrivals_agg');
matlabdatafile.veharrivals_CbyC = fullfile(datafolder,'veharrivals_CbyC');

% Load data from csv if not loaded already -----------------------------
if( ~exist([matlabdatafile.G '.mat'],'file') || ...
    ~exist([matlabdatafile.veharrivals_agg '.mat'],'file') || ...
    ~exist([matlabdatafile.veharrivals_CbyC '.mat'],'file') )

    numIntersections = 4;
    directions = {'NB','SB'};
    
    for i=1:numIntersections
        for j=1:length(directions)
            
            [i j]
            
            d = directions{j};
            name = [num2str(i) d];
            
            % read the data files .........................................
            N = xlsread([datafolder filesep name '.csv']);
            ind = find(isnan(N(:,1)));
            sectionindex = ind(diff(ind)~=1);
            clear ind
            
            % section 1: *** Probability of Green - Aggregated Distributions ***
            data = N(sectionindex(1)+1:sectionindex(2)-3,:);
            G(i).(d).time = data(:,1);
            G(i).(d).data = data(:,3:end);
            
            % section 2: *** Vehicle Arrivals - Aggregated Distributions ***'
            data = N(sectionindex(2)+1:sectionindex(3)-3,:);
            veharrivals_agg(i).(d).time = data(:,1);
            veharrivals_agg(i).(d).data = data(:,3:end);
            
            % section 3: *** Vehicle Arrivals - Cycle by Cycle ***'
            data = N(sectionindex(3)+1:end,:);
            veharrivals_CbyC(i).(d).time = data(:,1);
            veharrivals_CbyC(i).(d).data = data(:,3:end);
            
            clear N data
        end
    end
    
    save(matlabdatafile.G,'G')
    save(matlabdatafile.veharrivals_agg,'veharrivals_agg')
    save(matlabdatafile.veharrivals_CbyC,'veharrivals_CbyC')

end

% xxx -----------------------------

load(matlabdatafile.G)

figure
plot([1:114],G(1).NB.data(14,:))
set(gca,'XLim',[0 114])
grid
title(num2str(i))



disp('done')
