%% CoSMoS-AK: water level validation
% v1.0  Nederhoff   Jun-19
% v1.1  Nederhoff   May-20
% v1.2  Nederhof    2020-06-01
% v1.3  Nederhoff   2020-06-22      applicable with more years
clear all
close all
clc

% Get direction
maindir = 'q:\Projects\Alaska\CoMoS_AK\';

% Load observations (high-frequency and low-frequency)
cd(maindir)
load('01_data\waterlevels\coops\observations_refined.mat');
meter_feet          = 3.28084;

% Settings
spinup_time         = 7;                    % in days

% years wanted
years_wanted        = [2011:2018];

%% Several models
clear his_nc
his_nc_name{1}       = ['version006\normalruns\'];

for uniques = 1:length(his_nc_name)
    
    % Get version for this one
    version             = his_nc_name{uniques};
    
    %% 1. Retrieve data: takes 2 hours
    tic
    clear data data2
    cd(maindir)
    his_nc          = ['04_modelruns\', version, '\results_cluster\', num2str(years_wanted(1)), '\cosmos_ak_his.nc'];
    station_name    = nc_varget(his_nc, 'station_name');

    for ii = 1:length(observations)
        
        % Change name
        disp(['working on ' num2str(ii), ' of ', num2str(length(observations))]);
        for jj = 1:length(station_name)
            try
                name_TMP    = station_name(jj,:);
                name_TMP    = name_TMP(1:7);
                name_TMP    = name_TMP(find(~isspace(name_TMP)));
                name_TMP    = str2num(name_TMP);
                idfind(jj)  = str2num(observations(ii).IDcode) == name_TMP;
            catch
            end
        end
        idfind          = find(idfind == 1);
        
        if ~isempty(idfind)
            
            % Save data
            data(ii).IDcode             = observations(ii).IDcode;
            data(ii).Name               = observations(ii).Name;
            data(ii).lat                = observations(ii).y;
            data(ii).model.datetime     = [];
            data(ii).model.waterlevel   = [];
            data(ii).model.x            = nc_varget(his_nc, 'station_x_coordinate', [idfind-1], [1]);
            data(ii).model.y            = nc_varget(his_nc, 'station_y_coordinate', [idfind-1], [1]);
            
            % Get data of other years
            for yr = 1:length(years_wanted)
                
                % Define new his file
                clear xTMP yTMP
                his_nc                      = ['04_modelruns\', version, '\results_cluster\', num2str(years_wanted(yr)), '\cosmos_ak_his.nc'];
                reference_time              = ncreadatt(his_nc,'time','units');
                idtime                      = strfind(reference_time, 'since ');
                reference_time              = reference_time(:,[idtime+6: idtime+15]);
                reference_time              = datenum(reference_time, 'yyyy-mm-dd');
                xTMP                        = nc_varget(his_nc, 'time')/3600/24 + reference_time;
                yTMP                        = nc_varget(his_nc, 'waterlevel', [0 idfind-1], [Inf 1]);
                idwanted_model              = xTMP >= (min(xTMP) + spinup_time) & xTMP <= max(xTMP);

                % Model
                data(ii).model.datetime     = [data(ii).model.datetime  xTMP(idwanted_model)'];
                data(ii).model.waterlevel   = [data(ii).model.waterlevel  yTMP(idwanted_model)'];
            end
       
            % Retrieve data
            idwanted_data               = observations(ii).datetime >= min(data(ii).model.datetime) & observations(ii).datetime <= max(data(ii).model.datetime);
            data(ii).obs.datetime       = observations(ii).datetime(idwanted_data);
            data(ii).obs.waterlevel     = observations(ii).waterlevel(idwanted_data);
            data(ii).obs.coef_full      = observations(ii).coef;

            % Interpolate model results to same datetime as observed
            data(ii).model.waterlevel   = interpolation_wavetimeseries(data(ii).model.datetime, data(ii).model.waterlevel, data(ii).obs.datetime, (2/24));
            data(ii).model.datetime     = data(ii).obs.datetime;
            
            % Skill score
            [stat]                      = model_skill2(data(ii).obs.waterlevel, data(ii).model.waterlevel, data(ii).model.waterlevel);
            data(ii).model.skill        = stat;
            
            % Tide analysis
            % Model
            if ~isempty(data(ii).model.datetime) || ~isnan(nanmean(data(ii).model.waterlevel))
                
                % Determine tide and NTR on small subset
                coef                        = ut_solv(data(ii).model.datetime, data(ii).model.waterlevel, [],  data(ii).lat , 'auto', 'LinCI', 'OLS', 'White');
                data(ii).model.coef         = coef;
                data(ii).model.coef.mean    = 0;
                data(ii).model.tide         = ut_reconstr(data(ii).obs.datetime, data(ii).model.coef)';
                data(ii).model.tide         = data(ii).model.tide;
                for xx = 1:3
                    if xx == 1
                        ntr                                     = data(ii).model.waterlevel -  data(ii).model.tide'; 
                        jdf                                     = data(ii).model.datetime;
                    end
                    if xx == 2 & length(data(ii).model.datetime) > 2000
                        [ntr, jdf]                              = cmglowpass(data(ii).model.waterlevel, 360, data(ii).model.datetime,(1/24),'osu');
                    end
                    if xx == 3
                        dtwanted_1                              = (1/24/6); dtwanted_2 = 1/dtwanted_1;
                        xTMP                                    = nanmin(data(ii).model.datetime):dtwanted_1:nanmax(data(ii).model.datetime);
                        yTMP                                    = interp1(data(ii).model.datetime, data(ii).model.waterlevel, xTMP);
                        yTMP_smooth                             = smoothdata(yTMP, 'movmean', dtwanted_2, 'includenan');
                        ntr                                     = yTMP_smooth;
                        jdf                                     = xTMP;
                    end
                    data(ii).model.NTR(xx).values               = ntr;
                    data(ii).model.NTR(xx).datetime             = jdf;
                end
                
                % Data
                coef                        = ut_solv(data(ii).obs.datetime, data(ii).obs.waterlevel, [],  data(ii).lat , 'auto', 'LinCI', 'OLS', 'White');
                data(ii).obs.coef           = coef;
                data(ii).obs.coef.mean      = 0;
                data(ii).obs.tide           = ut_reconstr(data(ii).obs.datetime, data(ii).obs.coef)';
                data(ii).obs.tide           = data(ii).obs.tide;
                for xx = 1:3
                    if xx == 1
                        ntr                             = data(ii).obs.waterlevel -  data(ii).obs.tide'; 
                        jdf                             = data(ii).obs.datetime;
                    end
                    if xx == 2 & length(data(ii).obs.datetime) > 2000
                        [ntr, jdf]                      = cmglowpass(data(ii).obs.waterlevel, 360, data(ii).obs.datetime,(1/24),'osu');
                    end
                    if xx == 3
                        dtwanted_1              = (1/24/6); dtwanted_2 = 1/dtwanted_1;
                        xTMP                    = nanmin(data(ii).obs.datetime):dtwanted_1:nanmax(data(ii).obs.datetime);
                        yTMP                    = interp1(data(ii).obs.datetime, data(ii).obs.waterlevel, xTMP);
                        yTMP_smooth             = smoothdata(yTMP, 'movmean', dtwanted_2, 'includenan');    % daily averaged
                        ntr                     = yTMP_smooth;
                        jdf                     = xTMP;
                    end
                    data(ii).obs.NTR(xx).values                 = ntr;
                    data(ii).obs.NTR(xx).datetime               = jdf;
                end
                
                % Compute skill of tide and NTR
                [stat]                      = model_skill2(data(ii).obs.tide, data(ii).model.tide, data(ii).model.tide);
                data(ii).model.skill_tide   = stat;
                for xx = 1:3
                    [stat]                      = model_skill2(data(ii).obs.NTR(xx).values, data(ii).model.NTR(xx).values, data(ii).obs.NTR(xx).values);
                    data(ii).model.skill_NTR(xx)= stat;
                end
            end
        end
    end
    toc

    % Only continue with non-nan
    clear idnotnan
    for ii = 1:length(data)
        % Do we have a nan?
        try
            idnotnan(ii) = data(ii).model.skill.rmse;
        catch
            idnotnan(ii) = NaN;
        end
    end
    idnotnan    = ~isnan(idnotnan);
    data2       = data(idnotnan);
    data        = data2;
    
    %% 2. Plotting
    Cliner = linspecer(6);
    
    % Make space for amp/phi
    close all
    amps    = zeros(length(data),2,10);
    phis    = zeros(length(data),2,10);
    
    for ii = 1:length(data)
        
        % Go to folder
        cd(maindir)
        cd('05_post_processing\');
        mkdir(version); cd(version);
        mkdir('coops_wl'); cd('coops_wl');
        
        try
            
            %%  Tidal comparison
            close all
            Y = 29.7/3;   X = 21.0;
            xSize = X - 2*0.5;   ySize = Y - 2*0.5; % figure size on paper (width & height)
            hFig = figure('Visible','Off');
            hold on;
            set(hFig, 'PaperUnits','centimeters');
            set(hFig, 'PaperSize',[X Y]);
            set(hFig, 'PaperPosition',[0.5 0.5 xSize ySize]);
            set(hFig, 'PaperOrientation','portrait');
            
            % Components
            components_wanted = data(ii).obs.coef_full.name;
            components_wanted = components_wanted(1:10);
            for jj = 1:length(components_wanted)
                
                % Find matching onces for model
                idwanted                  = strmatch(components_wanted{jj},data(ii).model.coef.name,'exact');
                if ~isempty(idwanted)
                    amps(ii, 2,jj)  = data(ii).model.coef.A(idwanted);
                    ampslow(ii,2,jj)= data(ii).model.coef.A_ci(idwanted);
                    phis(ii,2,jj)   = data(ii).model.coef.g(idwanted);
                    phislow(ii,2,jj)= data(ii).model.coef.g_ci(idwanted);
                end
                
                % Find matching onces for observation
                idwanted                  = strmatch(components_wanted{jj},data(ii).obs.coef.name,'exact');
                if ~isempty(idwanted)
                    amps(ii, 1,jj)  = data(ii).obs.coef.A(idwanted);
                    ampslow(ii,1,jj)= data(ii).obs.coef.A_ci(idwanted);
                    phis(ii,1,jj)   = data(ii).obs.coef.g(idwanted);
                    phislow(ii,1,jj)= data(ii).obs.coef.g_ci(idwanted);
                end
            end
            amp_TMP = squeeze(amps(ii,:,:));
            phi_TMP = squeeze(phis(ii,:,:));
            amp_err = squeeze(ampslow(ii,:,:));
            phi_err = squeeze(phislow(ii,:,:));

            sub1 = subplot(2,1,1); hold on;
            hbar = bar(amp_TMP'); title([data(ii).Name, ' - ' , data(ii).IDcode]);
            set(hbar(1), 'FaceColor', Cliner(1,:), 'edgecolor','none')
            set(hbar(2), 'FaceColor', Cliner(2,:), 'edgecolor','none')
            set(sub1, 'xticklabel', '');
            box on; ylabel('tidal amplitude [meter]')

            % Get errors
            x = [];
            for i = 1:2
                x = [x ; hbar(i).XEndPoints];
            end
            errorbar(x',amp_TMP',(amp_err*1)','k','linestyle','none')
            
            legend('observed', 'modelled', '95% confidence intervals', 'orientation', 'horizontal')
            sub1.YLim = sub1.YLim*1.1;

            % Get phase
            yyaxis right
            ylabel(['tidal amplitude [feet]'])
            sub1.YAxis(2).Limits        = sub1.YAxis(1).Limits * meter_feet;
            sub1.YAxis(2).Color         = 'k';
            grid on; box on;
            sub2 = subplot(2,1,2); hold on;
            hbar = bar(phi_TMP');
            set(hbar(1), 'FaceColor', Cliner(1,:), 'edgecolor','none')
            set(hbar(2), 'FaceColor', Cliner(2,:), 'edgecolor','none')
            set(sub2, 'xticklabel', components_wanted);
            box on; ylabel('tidal amplitude [degrees]')
            grid on; box on;
            
            % Get errors
            x = [];
            for i = 1:2
                x = [x ; hbar(i).XEndPoints];
            end
            errorbar(x',phi_TMP',(phi_err*1)','k','linestyle','none')
            ylim([0 360])
           
            figure_png = ['tide', data(ii).IDcode, '.png'];
            print('-dpng','-r300',figure_png); close all
            
            %% Time series
            for ijn = 1:2
                
                % Create frame
                if ijn == 1
                    xwanted     = round(median(data(ii).obs.datetime));
                    xwanted     = [xwanted xwanted+7];
                    figure_png  = ['timeseries_random_', data(ii).IDcode, '.png'];
                    
                elseif ijn == 2
                    idmax       = find(max(data(ii).obs.waterlevel) == data(ii).obs.waterlevel); idmax = idmax(1);
                    idmax2      = find(data(ii).obs.datetime >= data(ii).obs.datetime(idmax)-3.5  & data(ii).obs.datetime <= data(ii).obs.datetime(idmax)+3.5);
                    xwanted     = [data(ii).obs.datetime(idmax)-3.5 data(ii).obs.datetime(idmax)+3.5];
                    figure_png  = ['timeseries_storm_', data(ii).IDcode, '.png'];
                end
                
                close all
                Y = 29.7/3;   X = 21.0;
                xSize = X - 2*0.5;   ySize = Y - 2*0.5; % figure size on paper (width & height)
                hFig = figure('Visible','Off');
                hold on;
                set(hFig, 'PaperUnits','centimeters');
                set(hFig, 'PaperSize',[X Y]);
                set(hFig, 'PaperPosition',[0.5 0.5 xSize ySize]);
                set(hFig, 'PaperOrientation','portrait');
                sub1 = subplot(1,1,1); hold on;
                
                % Show time series
                clear hplot
                hplot(3) = plot(data(ii).obs.datetime, data(ii).obs.tide + nanmean(data(ii).model.waterlevel), '--');
                hplot(2) = plot(data(ii).model.datetime, data(ii).model.waterlevel);
                hplot(1) = plot(data(ii).obs.datetime, data(ii).obs.waterlevel, '-.');
                %hplot(3) = plot(data(ii).obs.datetime, data(ii).obs.waterlevel-data(ii).model.waterlevel);
                
                set(hplot(1), 'Color', Cliner(1,:), 'linewidth', 1.5);
                set(hplot(2), 'Color', Cliner(2,:), 'linewidth', 1.5);
                set(hplot(3), 'Color', Cliner(4,:), 'linewidth', 1.0);
                title([data(ii).Name, ' - ' , data(ii).IDcode]);
                
                % Feet axis
                xlim([xwanted]);
                sub1.YLim = sub1.YLim*1.1;
                ylabel('water level (z_s) [meter + MSL')
                yyaxis right
                ylabel(['water level (z_s) [feet + MSL]'])
                
                % Time frame
                box on;
                hlegend  = legend([hplot], 'observed', 'modelled', 'tide', 'orientation', 'horizontal', 'location', 'southoutside');
                datetick('x','yyyy/mm/dd', 'keeplimits');

                if diff(sub1.YAxis(1).Limits) < 10
                    sub1.YAxis(2).Limits        = sub1.YAxis(1).Limits * meter_feet;
                    sub1.YAxis(2).Color         = 'k';
                    sub1.YAxis(2).TickValues    = round(sub1.YAxis(1).Limits(1) * meter_feet):1:round(sub1.YAxis(1).Limits(2) * meter_feet);
                end
                grid on; box on;
                print('-dpng','-r300',figure_png);
                close all
            end
%             
%             %% QQ-plot
%             addpath('c:\Matlab\personal\IoSR-Surrey-MatlabToolbox-4bff1bb\')
%             iosr.install
%             
%             close all
%             figure; hold on;
%             idwanted = ~isnan(data(ii).obs.waterlevel) & ~isnan(data(ii).model.waterlevel);
%             iosr.statistics.qqPlot([data(ii).obs.waterlevel(idwanted); data(ii).model.waterlevel(idwanted)]','both')
            
            %% Plot with fading
            Y = 29.7/2;   X = 21.0;
            xSize = X - 2*0.5;   ySize = Y - 2*0.5; % figure size on paper (width & height)
            hFig = figure('Visible','Off');
            hold on;
            set(hFig, 'PaperUnits','centimeters');
            set(hFig, 'PaperSize',[X Y]);
            set(hFig, 'PaperPosition',[0.5 0.5 xSize ySize]);
            set(hFig, 'PaperOrientation','portrait');
            hplot = scatter(data(ii).obs.waterlevel*1, data(ii).model.waterlevel*1, 'filled');
            set(hplot, 'Cdata', Cliner(3,:));    set(hplot, 'MarkerFaceAlpha', 0.05);
            hplot = plot([min(data(ii).obs.waterlevel)*1 max(data(ii).obs.waterlevel)*1], [min(data(ii).obs.waterlevel)*1 max(data(ii).obs.waterlevel)*1], '-k');
            xlim([min(data(ii).obs.waterlevel)*1 max(data(ii).obs.waterlevel)*1]);
            ylim([min(data(ii).obs.waterlevel)*1 max(data(ii).obs.waterlevel)*1])
            title([data(ii).Name, ' - ' , data(ii).IDcode]);
            xlabel('water level observed [m + MSL]');
            ylabel('water level measured [m + MSL]');
            htext1 = ['MAE_u: ', num2str( round(data(ii).model.skill.umae*1*100)/100), ' meter'];
            htext2 = ['bias: ', num2str( round(data(ii).model.skill.bias*1*100)/100), ' meter'];
            htext3 = ['RMSE: ', num2str( round(data(ii).model.skill.rmse*1*100)/100), ' meter'];
            text1  = text(0.1, 0.80, htext1, 'sc');
            text2  = text(0.1, 0.85, htext2, 'sc');
            text3  = text(0.1, 0.90, htext3, 'sc');
            box on;
            figure_png = ['scatter', data(ii).IDcode, '.png'];
            print('-dpng','-r300',figure_png); close all
            
            
            %% Error and NTR over time (more for me)
            close all
            Y = 29.7/3;   X = 21.0;
            xSize = X - 2*0.5;   ySize = Y - 2*0.5; % figure size on paper (width & height)
            hFig = figure%('Visible','Off');
            hold on;
            set(hFig, 'PaperUnits','centimeters');
            set(hFig, 'PaperSize',[X Y]);
            set(hFig, 'PaperPosition',[0.5 0.5 xSize ySize]);
            set(hFig, 'PaperOrientation','portrait');
            
            % Show time series
            %hplot = plot(boundary.in.time, boundary.in.wl*-1, '-k');
            error1 = data(ii).obs.waterlevel-data(ii).model.waterlevel;
            if median(diff(data(ii).obs.datetime)*24) < 0.5
                error2 = smoothdata(error1,'movmean', round(1/median(diff(data(ii).obs.datetime))), 'includenan');
            else
                error2 = smoothdata(error1,'movmean', round(1/median(diff(data(ii).obs.datetime))), 'includenan');
            end
            minerror = min(error1);
            maxerror = max(error1);
            range    = max([minerror*-1, maxerror]);
            
            hplot(4) = scatter(data(ii).obs.datetime, error1, 'filled');
            ylim([nanmin(error2)-0.2 nanmax(error2)+0.2]);
            set(hplot(4),'MarkerFaceAlpha', 0.1)
            set(hplot(4), 'MarkerFaceColor', Cliner(4,:));
            
            hplot(3) = plot(data(ii).obs.datetime, error2); hold on;
            set(hplot(3), 'Color', Cliner(3,:), 'linewidth', 2.0);
            
            
            hplot6 = plot([min(data(ii).obs.datetime) max(data(ii).obs.datetime)], [0 0], '--r');
            
            grid on; box on;
            title([data(ii).Name, ' - ' , data(ii).IDcode]);
            axis tight;
            datetick('x', 'yyyy-mm-dd', 'keeplimits');
            legend([hplot(4) hplot(3)], 'instanteous', 'daily moving average');
            ylabel('error [m]');             ylim([range*-1 range]);
            text(0.1, 0.75, 'observed higher', 'sc');
            text(0.1, 0.25, 'observed lower', 'sc');
            grid on; box on;
            figure_png = ['error_', data(ii).IDcode, '.png'];
            print('-dpng','-r300',figure_png); close all
            
            %% Make time series of NTR
            for ijn = 1:3
                
                close all
                Y = 29.7/3;   X = 21.0;
                xSize = X - 2*0.5;   ySize = Y - 2*0.5; % figure size on paper (width & height)
                hFig = figure('Visible','Off');
                hold on;
                set(hFig, 'PaperUnits','centimeters');
                set(hFig, 'PaperSize',[X Y]);
                set(hFig, 'PaperPosition',[0.5 0.5 xSize ySize]);
                set(hFig, 'PaperOrientation','portrait');
                sub1 = subplot(1,1,1); hold on;
                
                % Show time series
                clear hplot
                hplot(2) = plot(data(ii).model.NTR(ijn).datetime, data(ii).model.NTR(ijn).values);
                hplot(1) = plot(data(ii).obs.NTR(ijn).datetime, data(ii).obs.NTR(ijn).values, '--');
                set(hplot(1), 'Color', Cliner(1,:), 'linewidth', 1.5);
                set(hplot(2), 'Color', Cliner(2,:), 'linewidth', 1.5);
                title([data(ii).Name, ' - ' , data(ii).IDcode]);
                
                % Feet axis
                ylabel('water level (z_s) [meter + MSL')
                yyaxis right
                ylabel(['water level (z_s) [feet + MSL]'])
                
                % Time frame
                box on;
                hlegend  = legend(hplot,'observed', 'modelled', 'orientation', 'horizontal', 'location', 'southoutside');
                datetick('x','yyyy/mm/dd', 'keeplimits');
                
                sub1.YAxis(2).Limits        = sub1.YAxis(1).Limits * meter_feet;
                sub1.YAxis(2).Color         = 'k';
                sub1.YAxis(2).TickValues    = round(sub1.YAxis(1).Limits(1) * meter_feet):0.25:round(sub1.YAxis(1).Limits(2) * meter_feet);
                grid on; box on;
                
                figure_png = ['NTR', num2str(ijn), '_', data(ii).IDcode, '.png'];
                print('-dpng','-r300',figure_png); close all
                close all
            end
            
            %% Make 1 overview plot
            close all;
            Y = 29.7/3;   X = 21.0;
            xSize = X - 2*0.5;   ySize = Y - 2*0.5; % figure size on paper (width & height)
            hFig = figure('Visible','Off'); hold on;
            set(hFig, 'PaperUnits','centimeters');
            set(hFig, 'PaperSize',[X Y]);
            set(hFig, 'PaperPosition',[0.5 0.5 xSize ySize]);
            set(hFig, 'PaperOrientation','portrait');
            
            % Total
            fname_title = [data(ii).Name, ' - ' , data(ii).IDcode];
            sub1 = subplot(3,3,1); hold on;
            clear hplot
            hplot(2) = plot(data(ii).model.datetime, data(ii).model.waterlevel);
            hplot(1) = plot(data(ii).obs.datetime, data(ii).obs.waterlevel, ' . ');
            set(hplot(1), 'Color', Cliner(1,:), 'linewidth', 1.5);
            set(hplot(2), 'Color', Cliner(2,:), 'linewidth', 1.5);
            xwanted  = round(median(data(ii).obs.datetime));
            xwanted  = [xwanted xwanted+7];
            axis tight;  xlim([xwanted]);
            datetick('x', 'mmm/dd', 'keeplimits');
            ylabel('elevation [MSL+m]');
            box on; grid on;        title(fname_title);
            htext1 = ['RMSE: ' num2str(data(ii).model.skill.rmse*100,'%.1f'), ' [cm]'];
            htext2 = ['MAEu: ' num2str(data(ii).model.skill.umae*100,'%.1f'), ' [cm]'];
            htext3 = ['SCI: ' num2str(data(ii).model.skill.sci*100,'%.1f'), ' [%]'];
            text(0.1, 0.9, htext1, 'sc', 'Fontsize', 8);
            text(0.4, 0.9, htext2, 'sc', 'Fontsize', 8);
            text(0.7, 0.9, htext3, 'sc', 'Fontsize', 8);
            sub1.XAxis.TickValues = linspace(sub1.XAxis.Limits(1), sub1.XAxis.Limits(2), 8);
            for qaz = 1:length(sub1.XAxis.TickValues)
                wanted{qaz} = datestr(sub1.XAxis.TickValues(qaz), 'mmm/dd');
            end
            set(sub1, 'xticklabel', wanted)
            
            % Error scatter
            sub2 = subplot(3,3,2); hold on;
            hplot = scatter(data(ii).obs.waterlevel, data(ii).model.waterlevel, 'filled');   set(hplot, 'Cdata', Cliner(3,:));    set(hplot, 'MarkerFaceAlpha', 0.05);
            hplot = scatter(data(ii).obs.NTR(3).values, data(ii).model.NTR(3).values, 'filled');   set(hplot, 'Cdata', Cliner(4,:));    set(hplot, 'MarkerFaceAlpha', 0.05);
            hplot = plot([min(data(ii).obs.waterlevel) max(data(ii).obs.waterlevel)], [min(data(ii).obs.waterlevel) max(data(ii).obs.waterlevel)], '-k');
            xlabel('observed [m]');  ylabel('modelled [m]');
            grid on; box on;
            htext1 = ['R^2: ' num2str(data(ii).model.skill.r2,'%.2f'), ' [-]'];
            htext2 = ['R^2: ' num2str(data(ii).model.skill_NTR(2).r2,'%.2f'), ' [-]'];
            text(0.05, 0.85, htext1, 'sc', 'color', Cliner(3,:), 'Fontsize', 8);
            text(0.60, 0.1, htext2, 'sc', 'color', Cliner(4,:), 'Fontsize', 8);
            
            % Filtered
            sub3 = subplot(3,3,3); hold on;
            hplot(2) = plot(data(ii).model.NTR(3).datetime, data(ii).model.NTR(3).values);
            hplot(1) = plot(data(ii).obs.NTR(3).datetime, data(ii).obs.NTR(3).values, '--');
            set(hplot(1), 'Color', Cliner(1,:), 'linewidth', 1.5);
            set(hplot(2), 'Color', Cliner(2,:), 'linewidth', 1.5);
            datetick('x', 'mmm/dd', 'keeplimits');
            ylabel('elevation [MSL+m]');
            box on; grid on; axis tight
            htext1 = ['RMSE: ' num2str(data(ii).model.skill_NTR(end).rmse*100,'%.1f'), ' [cm]'];
            htext2 = ['MAEu: ' num2str(data(ii).model.skill_NTR(end).umae*100,'%.1f'), ' [cm]'];
            htext3 = ['SCI: ' num2str(data(ii).model.skill_NTR(end).sci*100,'%.1f'), ' [%]'];
            text(0.1, 0.9, htext1, 'sc', 'Fontsize', 8);
            text(0.4, 0.9, htext2, 'sc', 'Fontsize', 8);
            text(0.7, 0.9, htext3, 'sc', 'Fontsize', 8);
            legend('modeled', 'observed',' location', 'best');
            sub3.XAxis.TickValues = linspace(sub3.XAxis.Limits(1), sub3.XAxis.Limits(2), 8);
            for qaz = 1:length(sub3.XAxis.TickValues)
                wanted{qaz} = datestr(sub3.XAxis.TickValues(qaz), 'yyyy/mm/dd');
            end
            set(sub3, 'xticklabel', wanted)
            
            yearswanted   = year(min(data(ii).model.NTR(2).datetime)):year(max(data(ii).model.NTR(2).datetime));
            for yr = 1:length(yearswanted)
                xwanted_loc(yr) = datenum(yearswanted(yr),1,1);
                xwanted_tick{yr}= datestr(xwanted_loc(yr), 'yyyy/mm/dd');
            end
            set(sub3, 'xtick', xwanted_loc, 'xticklabel', xwanted_tick)

            
            % Where is this?
            Gnetwork = dflowfm.readNet('q:\Projects\Alaska\CoMoS_AK\04_modelruns\version002\meshes\ak_Liv_net.nc');
            Gnetwork.face = [];
            sub4 = subplot(3,3,4); hold on;
            %surf(xcol, ycol, zcol, geo); shading flat
            %view(0, 90)
            hplot = dflowfm.plotNet(Gnetwork, 'node', {'.k','markersize',0.001}); hold on;
            hplot = plot(data(ii).model.x, data(ii).model.y, 'm*', 'markersize', 5);
            axis equal
            axis off
            
            % Tidal amplitude
            if data(ii).obs.coef.A(1) > 0.1
                
                sub5 = subplot(3,3,5);
                hbar = bar(amp_TMP');
                ylim([0 max(max(amp_TMP))*1.05]);
                set(hbar(1), 'FaceColor', Cliner(1,:));
                set(hbar(2), 'FaceColor', Cliner(2,:));
                box on; ylabel('amplitude [meter]')
                set(sub5, 'xticklabel', components_wanted);
                xlim([0.5 5.5])
                
                sub6 = subplot(3,3,6);
                hbar = bar(phi_TMP');
                ylim([0 360]);
                set(hbar(1), 'FaceColor', Cliner(1,:));
                set(hbar(2), 'FaceColor', Cliner(2,:));
                box on; ylabel('phase [\circ]')
                set(sub6, 'xticklabel', components_wanted);
                xlim([0.5 5.5]); set(sub6, 'ytick', [0 90 180 270 360]);
            end
            
            % Set rest correct
            set(sub1, 'Position', [0.1 0.7 0.6 0.2]);       % time series of water level
            set(sub4, 'Position', [0.7 0.7 0.2 0.2])        % location
            
            set(sub2, 'Position', [0.1 0.4 0.2 0.2])        % scatter
            set(sub3, 'Position', [0.1 0.1 0.8 0.2])        % time series of tidally-filtered water level
            
            if data(ii).obs.coef.A(1) > 0.1
                set(sub5, 'Position', [0.4 0.4 0.2 0.2])        % amplitude
                set(sub6, 'Position', [0.7 0.4 0.2 0.2])        % ohase
            end
            figure_png = ['overview_', data(ii).IDcode, '.png'];
            print('-dpng','-r300',figure_png); close all
            
        catch
            disp(['something went wrong with ' data(ii).IDcode])
        end
    end
    
    %         %% 3. Overview plot
    %         % Get bathymetry
    %         % Get geoimage
    %         cd(maindir);
    %         [xcol, ycol, zcol, col] = georeference_image('01_data\geoimage_SFBD_UTM10');
    %         changexy    = 1/1000;
    %         areainterest.x = [475 625];
    %         areainterest.y = [4125 4250];
    %
    %         % Figure
    %         close all
    %         Y = 29.7/2;   X = 21.0;
    %         xSize = X - 2*0.5;   ySize = Y - 2*0.5; % figure size on paper (width & height)
    %         hFig = figure('Visible','Off');
    %         hold on;
    %         set(hFig, 'PaperUnits','centimeters');
    %         set(hFig, 'PaperSize',[X Y]);
    %         set(hFig, 'PaperPosition',[0.5 0.5 xSize ySize]);
    %         set(hFig, 'PaperOrientation','portrait');
    %         surf(xcol*changexy, ycol*changexy, zcol, col); shading flat;
    %         xlabel('x in UTM 10N [km]'); ylabel('y in UTM 10N [km]');
    %         grid on; box on; xlim(areainterest.x); ylim(areainterest.y)
    %         axis equal
    %         xlim(areainterest.x); ylim(areainterest.y)
    %
    %         % Plotting
    %         for jj = 1:length(data)
    %
    %             % Location
    %             hcircle(jj) = plot(data(jj).model.x*changexy, data(jj).model.y*changexy, 'wo');
    %             htext(jj) = text(data(jj).model.x*changexy+1, data(jj).model.y*changexy, data(jj).Name); set(htext(jj), 'color', 'w');
    %
    %         end
    %         figure_png = ['overview.png'];
    %         print('-dpng','-r300',figure_png); close all
    
    %% 4. Skill scores
    clear writing
    for jj = 1:length(data)
        writing{jj+1,1} = data(jj).IDcode;
        writing{jj+1,2} = data(jj).Name;
        
        try
            writing{jj+1,3} = num2str(data(jj).model.skill.rmse*100, '%2.1f');
            writing{jj+1,4} = num2str(data(jj).model.skill.urmse*100, '%2.1f');
            writing{jj+1,6} = num2str(data(jj).model.skill_NTR(1).rmse*100, '%2.1f');
            writing{jj+1,7} = num2str(data(jj).model.skill_NTR(2).rmse*100, '%2.1f');
            writing{jj+1,8} = num2str(data(jj).model.skill_HW.rmse*100, '%2.1f');
            writing{jj+1,9} = num2str(data(jj).model.skill_waterlevelMSL.rmse*100, '%2.1f');
            writing{jj+1,10} = num2str(data(jj).model.skill.mae*100, '%2.1f');
            writing{jj+1,11} = num2str(data(jj).model.skill.bias*100, '%2.1f');
        catch
        end
        
        try
            writing{jj+1,5} = num2str(data(jj).model.skill_tide.rmse*100, '%2.1f');
        catch
        end
        
    end
    writing{1,1} = 'IDcode';
    writing{1,2} = 'Name';
    writing{1,3} = 'RMSE [cm]';
    writing{1,4} = 'RMSE_u [cm]';
    writing{1,5} = 'RMSE tide [cm]';
    writing{1,6} = 'RMSE NTR (1 - Kees) [cm]';
    writing{1,7} = 'RMSE NTR (2 - Babak) [cm]';
    writing{1,8} = 'RMSE HW [cm]';
    writing{1,9} = 'RMSE WL [cm]';
    writing{1,10} = 'MAE [cm]';
    writing{1,11} = 'bias [cm]';
    
    cd(maindir)
    cd('05_post_processing\');
    mkdir(version); cd(version);
    mkdir('coops_wl'); cd('coops_wl');
    delete('skill.xls')
    try
        xlswrite('skill.xls', writing);
    catch
        pause(10);
        xlswrite('skill.xls', writing);
    end
    save('results.mat', 'data');
    
    
    % Make KML with names
    for xx = 1:length(data)
        lon(xx)     = data(xx).model.x;
        lat(xx)     = data(xx).model.y;
        name{xx}    = [data(xx).Name ' - ', data(xx).IDcode];
    end
    KMLtext(lat,lon,name,'fileName','locations.kml');

end
disp('done!')
