%% USGS CoSMoS-AK: compare coops based observations compared to FES
% v1.0  Nederhoff   Jan-19
% v1.1  Nederhoff   Jun-19
% v1.2  Nederhoff   2020-06-01
clear all
close all
clc

%% Load data
load('q:\Projects\Alaska\CoMoS_AK\01_data\waterlevels\coops\observations_refined.mat');
destout     = 'q:\Projects\Alaska\CoMoS_AK\02_data_analysis\waterlevels\FES\';

%% 4. Restructuring data and do basic tidal analysis
% Not needed anymore

%% 5. Readme + figure of the water levels and compare to TPXO
components_wanted = {'SA';'SSA';'M2';'O1';'S2';'K1'};       % 10 largest based on long record
for ii = 1:length(components_wanted)
    
    % Read FES M2
    fnc                     = ['c:\DUSA\_other\data\FES\', components_wanted{ii}, '.nc'];
    lon                     = nc_varget(fnc, 'lon');
    lat                     = nc_varget(fnc, 'lat');
    [lon, lat]              = meshgrid(lon,lat);
    tidal_amplitude_h       = nc_varget(fnc, 'amplitude');
    tidal_phase_h       	= nc_varget(fnc, 'phase');

    % Fix
    tidal_amplitude_h       = tidal_amplitude_h/100;
    rad                     = deg2rad(tidal_phase_h);
    idflip                  = rad < 0;
    rad(idflip)             = rad(idflip) + 2*pi;
    tidal_phase_h           = double(rad2deg(rad));
        
    % Figure amplitude
    close all
    A4fig
    hold on;
    hpcolor = pcolor(lon-360, lat, tidal_amplitude_h); shading flat
    axis equal
    xlim([-205 -134])
    ylim([48 81]);
    clear hscatter htext
    datasaved = NaN(length(observations),1);
    for jj = 1:length(observations)
        try
        idfindsave = [];
        for nn = 1:length(observations(jj).coef.name)
            idfind = strcmpi(observations(jj).coef.name{nn}, components_wanted{ii});
            if idfind == 1
                idfindsave = nn;
            end
        end
        hscatter(jj)    = scatter(observations(jj).x, observations(jj).y, [], observations(jj).coef.A(idfindsave), 'filled');
        %htext(jj)       = text(observations(jj).x+0.02, observations(jj).y, observations(jj).Name);
        %set(htext(jj), 'Color', 'w')
        set(hscatter(jj), 'MarkerEdgeColor' ,'w');
        datasaved(jj)   = hscatter(jj).CData;
        catch
        end
    end
    caxis([0 nanmax(datasaved)])
    grid on; box on;
    xlabel('longitude [\circ]')
    ylabel('latitude [\circ]')
    hcolorbar =colorbar;
    legend('FES', 'observed');
    ylabel(hcolorbar,components_wanted{ii});
    colormap(jet); title('Tidal consituents observed (NOAA) and FES-based')
    cd('q:\Projects\Alaska\CoMoS_AK\02_data_analysis\waterlevels\FES\')
    fname = ['FES_', components_wanted{ii}, '_amplitude.png'];
    print('-dpng','-r300', fname);
    close all
    
    % Figure phase
    close all
    A4fig
    hold on;
    hpcolor = pcolor(lon-360, lat, tidal_phase_h); shading flat
    axis equal
    xlim([-205 -134])
    ylim([48 81]);
    clear hscatter htext
    datasaved = NaN(length(observations),1);
    for jj = 1:length(observations)
        try
        idfindsave = [];
        for nn = 1:length(observations(jj).coef.name)
            idfind = strcmpi(observations(jj).coef.name{nn}, components_wanted{ii});
            if idfind == 1
                idfindsave = nn;
            end
        end
        hscatter(jj)    = scatter(observations(jj).x, observations(jj).y, [], observations(jj).coef.g(idfindsave), 'filled');
        %htext(jj)       = text(observations(jj).x+0.02, observations(jj).y, observations(jj).Name);
        %set(htext(jj), 'Color', 'w')
        set(hscatter(jj), 'MarkerEdgeColor' ,'w');
        datasaved(jj)   = hscatter(jj).CData;
        catch
        end
    end
    cmap = cmocean('phase');
    caxis([0 360])
    grid on; box on;
    xlabel('longitude [\circ]')
    ylabel('latitude [\circ]')
    hcolorbar =colorbar;
    legend('FES', 'observed');
    ylabel(hcolorbar,components_wanted{ii});
    colormap(cmap); 
    title('Tidal consituents observed (NOAA) and FES-based')
    cd('q:\Projects\Alaska\CoMoS_AK\02_data_analysis\waterlevels\FES\')
    fname = ['FES_', components_wanted{ii}, '_phase.png'];
    print('-dpng','-r300', fname);
    close all
    
end


