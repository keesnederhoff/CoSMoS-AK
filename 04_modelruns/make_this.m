%% Matlab script that prepares a model run
clear all
close all
clc

% Settings
destin_MAIN = 'q:\Projects\Alaska\CoMoS_AK\03_modelsetup\version007\';
destout     = 'q:\Projects\Alaska\CoMoS_AK\04_modelruns\version007\normalruns\';

% Variations
WY          = [2011:2018];

%% Loop
countup = 0;
for ii =1:length(WY)
    
    % Get folder ready
    destout_TMP     = [destout, 'year', num2str(WY(ii))];
    try
        rmdir(destout_TMP, 's')
        disp('folder deleted');
    catch
    end
    mkdir(destout_TMP); cd(destout_TMP);
    copyfile(destin_MAIN,destout_TMP)
    
    % Times (full water years with 7 days of spin-up)
    WYnow       = WY(ii); disp(['  WY', num2str(WYnow)])
    reftime     = datenum(WYnow-1,12,1);
    spinup      = 7;
    time_start  = reftime-spinup;
    time_end    = datenum(WYnow+1,1,1);
    
    % Time sets
    dt_user     = 120;                          % dt_user in seconds
    dtfac       = 86400; 	                    % second in a day
    latitude    = 65;
        
    % Define meteo
    cd(destout_TMP);
    FORCINGWANTED   = ['ERA5_normal', num2str(WYnow), '.nc'];
    [succes]        = replace_text('forcing_old.ext','FORCINGWANTED', FORCINGWANTED);
        
    % Set times in mdu
    cd(destout_TMP);
    [succes]        = replace_text('cosmos_ak.mdu','STARTDATE', datestr(reftime, 'yyyymmdd'));
    [succes]        = replace_text('cosmos_ak.mdu','TSTART_WANTED', num2str((time_start - reftime)*dtfac) );
    [succes]        = replace_text('cosmos_ak.mdu','TEND_WANTED', num2str((time_end - reftime)*dtfac) );
    
    % Set fourier
    cd(destout_TMP);
    components  ={'M2','S2','N2','K2','K1','O1','P1','Q1'};
    [succes]    = determine_fourier(reftime,time_start, time_end, spinup, latitude, dt_user, components, 'cosmos_ak.fou');
    fclose('all');
    
    % Save name in one long list
    countup                     = countup+1;
    dirs_simulation{countup}    = destout_TMP;
    
end