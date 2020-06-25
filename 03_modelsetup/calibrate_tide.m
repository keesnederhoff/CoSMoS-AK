%% Determine change values for CoMoS-AK
% v1.0  Nederhoff   2020-06-25
clear all
close all
clc

% Get input of where the tidal boundaries are
load('q:\Projects\Alaska\CoMoS_AK\03_modelsetup\version007\boundaries_FES2014\tidal.mat')

% Load data
load('q:\Projects\Alaska\CoMoS_AK\05_post_processing\version006\normalruns\coops_wl\results.mat')

% Make pairs and determine fit (future same for north)
% NB: i tried a couple of things regarding WHICH consituents to use
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
            idwanted(1)            = strmatch(pair(ii).consituents{cc}, data_subset.obs.coef.name,'exact');
            idwanted(2)            = strmatch(pair(ii).consituents{cc}, data_subset.model.coef.name,'exact');
            amp_factor(jj,cc)      = data_subset.obs.coef.A(idwanted(1)) ./ data_subset.model.coef.A(idwanted(2));
            
            % Find phase difference
            phase_difference(jj,cc)= data_subset.obs.coef.g(idwanted(1)) - data_subset.model.coef.g(idwanted(2));
            
        end
    end
    
    % NaN out large changes
    if ii == 1
        
        % Not too extreme
        amp_factor(amp_factor>1.5 | amp_factor < 0.5) = NaN;
        phase_difference(phase_difference>45 | phase_difference < -45) = NaN;
        
        % We simply apply the mean
        amp_factor_mean             = nanmean(amp_factor);
        phase_difference_mean       = nanmean(phase_difference);
        
    else
        amp_factor_mean             = amp_factor;
        phase_difference_mean       = phase_difference;
    end
    
    % Create bc correction files
    for jj = 1:length(pair(ii).boundaries)
        
        % Get BND
        cd('q:\Projects\Alaska\CoMoS_AK\03_modelsetup\version007\boundaries_FES2014\');
        bnd          = dflowfm_io_bnd('read', pair(ii).boundaries{jj});
    
        % Make tide
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
end

