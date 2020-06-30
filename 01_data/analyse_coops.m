%% USGS CoSMoS-AK: convert to raw to processed water levels
% v1.0  Nederhoff   Jan-19
% v1.1  Nederhoff   Jun-19
% v1.2  Nederhoff   2020-06-29
clear all
close all
clc

%% 0. Settings
% co-ops mat and where
destout     = 'q:\Projects\Alaska\CoMoS_AK\01_data\waterlevels\coops\figures\';
load 'q:\Projects\Alaska\CoMoS_AK\01_data\waterlevels\coops\co-ops_RAW_20200629.mat'

%% 4. Restructuring data and do basic tidal analysis
cd(destout); 
observations = observations_saved;
for ii = 1:length(observations)

        % Get basic values
        disp(['Working on: ', num2str(ii), ' of ', num2str(length(observations))]);
        clear idTMP values times lat
        values    = observations(ii).waterlevel;
        times     = observations(ii).datetime;
        lat       = observations(ii).y;
        dtdate    = median(diff(times));
        OUT_time  = min(times):dtdate:max(times);
        [OUT_hs, OUT_hs_old] = interpolation_wavetimeseries(times, values, OUT_time', (dtdate*2));

        % Save data
        observations(ii).waterlevel = OUT_hs;
        observations(ii).datetime   = OUT_time;
        
        % Reduce time period
        idwanted                     = (observations(ii).datetime > datenum(1900,1,1));         % basically, all the data we have
        
        % Apply u_tide (can take a while if there is a bunch of data)
        coef                        = ut_solv(OUT_time(idwanted), OUT_hs(idwanted), [], lat, 'auto', 'LinCI', 'OLS', 'White');
        [yout, ~ ]                  = ut_reconstr(OUT_time(idwanted), coef);
        
        % Check: makes sense?
        close all
        Y = 29.7/3;   X = 21.0;                       
        xSize = X - 2*0.5;   ySize = Y - 2*0.5; % figure size on paper (width & height)
        hFig = figure; hold on;
        set(hFig, 'PaperUnits','centimeters');
        set(hFig, 'PaperSize',[X Y]);
        set(hFig, 'PaperPosition',[0.5 0.5 xSize ySize]);
        set(hFig, 'PaperOrientation','portrait');; hold on;
        plot(observations(ii).datetime, yout);
        plot(observations(ii).datetime, observations(ii).waterlevel, '--r');
        grid on; box on;
        legend('tide', 'data');
        xlim([datenum(2019,1,1) datenum(2019,1,7)])
        datetick('x', 'dd/mmm/yy', 'keeplimits');
        fname = [observations(ii).IDcode, '.png'];
        print('-dpng','-r300',fname);

        % Do t_predict (for entire time period)
        observations(ii).tide       = yout;
        observations(ii).coef       = coef;
end
    
%% 5. Readme + figure of the water levels and compare to TPXO
% Read TPXO
fnc                     = 'c:\DUSA\_other\software\DelftDashboard\working_OET\data\tidemodels\tpxo80\tpxo80.nc';
lon                     = nc_varget(fnc, 'lon');
lat                     = nc_varget(fnc, 'lat');
[lon, lat]              = meshgrid(lon,lat);
tidal_amplitude_h       = permute(nc_varget(fnc, 'tidal_amplitude_h'),[3,2,1]);
tidal_phase_h       	= permute(nc_varget(fnc, 'tidal_phase_h'), [3,2,1]);
tidal_constituents      = nc_varget(fnc, 'tidal_constituents');
depth                   = nc_varget(fnc, 'depth')';
idnotwanted             = depth == 0;
tidal_amplitude_h(1,idnotwanted) = NaN;
tidal_phase_h(1,idnotwanted) = NaN;

% Figure
close all
A4fig
hold on;
hpcolor = pcolor(lon-360, lat, squeeze(tidal_amplitude_h(2,:,:))); shading flat
axis equal
xlim([-205 -134])
ylim([48 81]);
set(hpcolor, 'Facealpha', 0.75);
for ii = 1:length(observations)
    idfindsave = [];
    for jj = 1:length(observations(ii).coef.name)
        idfind = strcmpi(observations(ii).coef.name{jj}, 'M2');
        if idfind == 1
            idfindsave = jj;
        end
    end
    hscatter(ii)    = scatter(observations(ii).x, observations(ii).y, [], observations(ii).coef.A(idfindsave), 'filled');
    htext(ii)       = text(observations(ii).x+0.02, observations(ii).y, observations(ii).Name);
    set(htext(ii), 'Color', 'k')
    set(hscatter(ii), 'MarkerEdgeColor' ,'w');
end
caxis([0 1]);
grid on; box on;
xlabel('longitude [\circ]')
ylabel('latitude [\circ]')
hcolorbar =colorbar;
ylabel(hcolorbar, 'M2 amplitude [m]')
caxis([0 1.0]);
colormap(jet);
fname = 'tpxo_stations_m2.png';
print('-dpng','-r300', fname);
close all

%% 6. Save data
cd('q:\Projects\Alaska\CoMoS_AK\01_data\waterlevels\coops\');
save('observations_refined.mat', 'observations');