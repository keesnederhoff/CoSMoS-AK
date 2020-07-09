%% Determine change values for CoMoS-AK
% v1.0  Nederhoff   2020-06-25
% v1.1  Nederhoff   2020-06-30
clear all
close all
clc

% Get input of where the tidal boundaries are, load data and define if recalibration
load('q:\Projects\Alaska\CoMoS_AK\03_modelsetup\version007\boundaries_FES2014\tidal.mat')
load('q:\Projects\Alaska\CoMoS_AK\05_post_processing\version007\results_cluster\attempt1_20200625\coops_wl\results.mat')
recalibrate_file    = 'tidal_cor_version001_20200625.mat';       % leave empty is you do not want to re-calibrate
destout             = 'q:\Projects\Alaska\CoMoS_AK\03_modelsetup\version007\boundaries_FES2014\';
fname_version       = 'version002_20200702';

% Make pairs and determine fit (future same for north)
% NB: I tried a couple of things regarding WHICH consituents to use
% Pair 1
pair(1).boundaries  = {'tide_south.pli'; 'tide_southeast.pli'; 'tide_southwest.pli'};
pair(1).idwanted    = {'9455090';   '9457292';      '9459450';  '9459881';  '9461380'; '9461710'; '9462450'; '9462620'};       
                       % Seaward,   Kodiak Island,  Sand Point  % King Cove
                       % Adak Island, Atka, Nikloski, Unalaska
pair(1).consituents = {'J1','K1','K2','M2','MF','N2','NU2','O1','P1','Q1','S2','SA','SSA'};

% Pair 2
pair(2).boundaries  = {'tide_north.pli'; 'tide_northeast.pli'; 'tide_northwest.pli'};
pair(2).idwanted    = {'9497645'};       
                       % Only Prudhoe Bay, because Red Dog is influenced by
                       % friction
pair(2).consituents = {'K1','M2','MF','N2','O1','Q1','S2','SA','SSA'};


%% Get differences and factors
for ii = 1:length(pair)
    
    % Part A: determine value
    longlist = [];
    clear amp_factor phase_difference
    for jj = 1:length(pair(ii).idwanted)
        
        % Find results (make subset)
        clear idfind
        for nn = 1:length(data)
            name_TMP    = data(nn).IDcode;
            name_TMP    = str2num(name_TMP);
            idfind(nn)  = str2num(pair(ii).idwanted{jj}) == name_TMP;
        end
        idfind          = find(idfind == 1);
        data_subset     = data(idfind);
        longlist        = [longlist data_subset.obs.coef.name(1:10)'];    
    
        % Get consituents
        for cc = 1:length(pair(ii).consituents)
            
            % Find in factor for amplitude
            clear idwanted
            try
                idwanted(1)            = strmatch(pair(ii).consituents{cc}, data_subset.obs.coef.name,'exact');
                idwanted(2)            = strmatch(pair(ii).consituents{cc}, data_subset.model.coef.name,'exact');
                amp_factor(jj,cc)      = data_subset.obs.coef.A(idwanted(1)) ./ data_subset.model.coef.A(idwanted(2));
                phase_difference(jj,cc)= data_subset.obs.coef.g(idwanted(1)) - data_subset.model.coef.g(idwanted(2));
            catch
                amp_factor(jj,cc)      = NaN;
                phase_difference(jj,cc)= NaN;
            end
        end
    end
    
    % NaN out large changes
    if ii == 1 || ii == 2
        
        % Not too extreme
        amp_factor(amp_factor>1.5 | amp_factor < 0.5) = NaN;
        phase_difference(phase_difference>20 | phase_difference < -20) = NaN;
        
        % We simply apply the mean
        if size(amp_factor,1) > 1
            amp_factor_mean             = nanmean(amp_factor);
            phase_difference_mean       = nanmean(phase_difference);
        else
            amp_factor_mean             = amp_factor;
            phase_difference_mean       = phase_difference;
        end
    else
        amp_factor_mean             = amp_factor;
        phase_difference_mean       = phase_difference;
    end
    
    % Determine combination of current and previous cor files
    if ~isempty(recalibrate_file)
        
        % Get previous values
        cd(destout)
        recalibrate = load(recalibrate_file);
        
        % Combine
        amp_factor_mean2            = amp_factor_mean.*recalibrate.results(ii).amp_factor_mean;
        phase_difference_mean2      = phase_difference_mean+recalibrate.results(ii).phase_difference_mean;

        % Replace NaN and SSA/SA (long-term simulations needed)recalibrate.results(ii).amp_factor_mean
        amp_factor_mean2(isnan(amp_factor_mean2))               = recalibrate.results(ii).amp_factor_mean(isnan(amp_factor_mean2));
        phase_difference_mean2(isnan(phase_difference_mean2))   = recalibrate.results(ii).phase_difference_mean(isnan(phase_difference_mean2));

        % But back SSA/SA
        amp_factor_mean2(end-1:end)         = recalibrate.results(ii).amp_factor_mean(end-1:end);
        phase_difference_mean2(end-1:end)   = recalibrate.results(ii).phase_difference_mean(end-1:end);
               
        % Done
        amp_factor_mean                     = amp_factor_mean2;
        phase_difference_mean               = phase_difference_mean2;
    end
    
    
    % Create bc correction files
    for jj = 1:length(pair(ii).boundaries)
        
        % Get BND
        cd(destout);
        bnd          = dflowfm_io_bnd('read', pair(ii).boundaries{jj});
    
        % Make tide
        fname        = ['cor_' fname_version];
        mkdir(fname); cd(fname);
        fout = [pair(ii).boundaries{jj}(1:end-4), '_cor.bc'];
        fid=fopen(fout,'wt');

        for bb = 1:length(bnd.DATA)

            % Write out bc file with astronomical components
            % Write header
            fprintf(fid,'%s\n',['[forcing]']);
            fprintf(fid,'%s\n',['Name                            = ' bnd.Name{bb}]);
            fprintf(fid,'%s\n',['Function                        = astronomic-correction']);
            fprintf(fid,'%s\n',['Quantity                        = astronomic component']);
            fprintf(fid,'%s\n',['Unit                            = -']);
            fprintf(fid,'%s\n',['Quantity                        = waterlevelbnd amplitude']);
            fprintf(fid,'%s\n',['Unit                            = m']);
            fprintf(fid,'%s\n',['Quantity                        = waterlevelbnd phase']);
            fprintf(fid,'%s\n',['Unit                            = deg']);

            % Write values
            for kk = 1:length(pair(ii).consituents)

                % Write stuff
                xvalueTMP     = upper(pair(ii).consituents{kk});
                yvalueTMP     = amp_factor_mean(kk);
                zvalueTMP     = phase_difference_mean(kk);

                % Done
                fprintf(fid,'%s %8.6f %8.6f \n',xvalueTMP,yvalueTMP, zvalueTMP);
            end
            fprintf(fid,'%s\n',[]);

        end
        fclose('all');
    end
    
    % Save calibration values for reference (for re-calibration mainly)
    results(ii).boundaries                  = pair(ii).boundaries;
    results(ii).idwanted                    = pair(ii).idwanted;
    results(ii).consituents                 = pair(ii).consituents;
    results(ii).amp_factor_mean             = amp_factor_mean;
    results(ii).phase_difference_mean       = phase_difference_mean;

end

% Save data
cd(destout);
save(['tidal_cor_', fname_version, '.mat'], 'results');

