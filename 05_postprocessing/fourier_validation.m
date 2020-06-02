%% Alaska: compare Fourier file with observations
% v1.0  Nederhoff   2020-06-01
clear all
close all
clc

% Get direction
maindir = 'q:\Projects\Alaska\CoMoS_AK\';

% Obs
load('q:\Projects\Alaska\CoMoS_AK\01_data\waterlevels\coops\observations_refined.mat')
components_wanted = {'M2';'O1';'S2';'K1'};       % 5 expected most important once

% Make KML for GIS
for ii = 1:length(observations)
    lat(ii)     = observations(ii).y;
    lon(ii)     = observations(ii).x;
    IDcode{ii}  = observations(ii).IDcode;
end
cd('q:\Projects\Alaska\CoMoS_AK\03_modelsetup\version003_Liv_WGS84\KML\')
KMLtext(lat,lon,IDcode,'fileName', 'coops.kml')

%% Several models
clear his_nc
his_nc_name{1}       = ['version004\run01_serial\'];

for uniques = 1:length(his_nc_name)
    
    % Get version for this one
    version             = his_nc_name{uniques};
    
    %% 1. Retrieve data
    clear data data2
    cd(maindir)
    his_nc          = ['04_modelruns\', version, '\DFM_OUTPUT_cosmos_ak\cosmos_ak_fou.nc'];

    % Get information from t_tide
    tt          = t_getconsts;
    names       = tt.name;
    freqs       = tt.freq;
    for i=1:size(names,1)
        cnsts{i}=deblank(names(i,:));
    end
    for jj = 1:length(components_wanted)
        freqwanted(jj) = freqs(strmatch(components_wanted{jj},cnsts,'exact'));
    end
    periodwanted = 360./(1./freqwanted);
    
    % Get information from D-Flow FM
    for jj = 1:9
        try
            fm_freq_found(jj)       = nc_attget(his_nc, ['mesh2d_fourier00', num2str(jj),'_amp'], 'Frequency_degrees_per_hour');
        catch
            fm_freq_found(jj)       = NaN;
        end
    end
    clear idfind
    for jj = 1:length(periodwanted)
        [index,distance] = near(fm_freq_found,periodwanted(jj),1);
        if distance < 0.01
            idfind(jj) = index;
        else
            idfind(jj) = NaN;
        end
    end

    
    %% 2. Load per consituent (for both FES and D-Flow FM)
    for ii = 1:length(components_wanted)

        % A. Read D-Flow FM
        data.G                      = dflowfm.readNet(his_nc);
        TMP1                        = nc_varget(his_nc, ['mesh2d_fourier00', num2str(idfind(ii)), '_amp']); 
        TMP2                        = nc_varget(his_nc, ['mesh2d_fourier00', num2str(idfind(ii)), '_phs']);
        idreplace                   = TMP1 < 0.00;
        TMP1(idreplace)             = NaN;  TMP2(idreplace)         = NaN; 
        data.consituent(ii).amp     = TMP1;
        data.consituent(ii).phs     = TMP2;
        
        % Find min_depth
        fileinfo                    = nc_info(his_nc);
        for jj = 1:length(fileinfo.Dataset)
            idfindmin(jj) = contains(fileinfo.Dataset(jj).Name, 'min_depth');
        end
        idfindmin                   = find(idfindmin);
        TMP1                        = nc_varget(his_nc, fileinfo.Dataset(idfindmin).Name);
        data.minwaterdepth          = TMP1;
        
        % B. Read FES M2
        fnc                         = ['c:\DUSA\_other\data\FES\', components_wanted{ii}, '.nc'];
        lon                         = nc_varget(fnc, 'lon');
        lat                         = nc_varget(fnc, 'lat');
        [lon, lat]                  = meshgrid(lon,lat);
        lon                         = lon -360;
        tidal_amplitude_h           = nc_varget(fnc, 'amplitude');
        tidal_phase_h               = nc_varget(fnc, 'phase');

        % Fix
        tidal_amplitude_h           = tidal_amplitude_h/100;
        rad                         = deg2rad(tidal_phase_h);
        idflip                      = rad < 0;
        rad(idflip)                 = rad(idflip) + 2*pi;
        tidal_phase_h               = double(rad2deg(rad));

        % Save some points across the domain
        idwanted                    = find(lon > min(data.G.node.x) & lon < max(data.G.node.x) & lat > min(data.G.node.y) & lat < max(data.G.node.y) );
        idwanted                    = idwanted(1:1000:end);
        fes(ii).lon                 = lon(idwanted);
        fes(ii).lat                 = lat(idwanted);
        fes(ii).amp                 = tidal_amplitude_h(idwanted);
        fes(ii).phs                 = tidal_phase_h(idwanted);

    end
    
    %% 3. Plotting per consituent
    cd(maindir)
    cd('05_post_processing\');
    mkdir(version); cd(version);
    mkdir('fourier'); cd('fourier');
    
    for ii = 1:length(components_wanted)
        for xx = 1:2
            close all
            Y = 29.7/3;   X = 21.0;                       
            xSize = X - 2*0.5;   ySize = Y - 2*0.5; % figure size on paper (width & height)
            hFig = figure; hold on;
            set(hFig, 'PaperUnits','centimeters');
            set(hFig, 'PaperSize',[X Y]);
            set(hFig, 'PaperPosition',[0.5 0.5 xSize ySize]);
            set(hFig, 'PaperOrientation','portrait');    hold on;

            % 2D of D-Flow FM
            G = data.G; 
            if xx== 1; D.face = data.consituent(ii).amp; end
            if xx== 2
                TMP                         = data.consituent(ii).phs;   
                rad                         = deg2rad(TMP);
                idflip                      = rad < 0;
                rad(idflip)                 = rad(idflip) + 2*pi;
                D.face                      = double(rad2deg(rad));
            end
            face.mask   = data.minwaterdepth > 1;
            tri.mask    = face.mask(G.map3);
            rangenumbers= sort(D.face(D.face>0));
            numbers     = quantille_numbers(rangenumbers);
            h           = trisurfcorcen(G.tri(tri.mask,:),G.node.x,G.node.y,D.face(G.map3(tri.mask)));
            set(h,'edgeColor','none'); 
            view(0,90);

            % Add FES2014 dots
            if xx== 1; hscatter = scatter(fes(ii).lon, fes(ii).lat, 10, fes(ii).amp, 'filled'); end
            if xx== 2; hscatter = scatter(fes(ii).lon, fes(ii).lat, 10, fes(ii).phs, 'filled'); end

            % Add observed stations (bigger circles)
            clear lat lon amp
            for jj = 1:length(observations)
                try
                    lat(jj)     = observations(jj).y;
                    lon(jj)     = observations(jj).x;
                    ID          = strmatch(components_wanted(ii), observations(jj).coef.name, 'exact');
                    amp(jj)     = observations(jj).coef.A(ID);
                    phs(jj)     = observations(jj).coef.g(ID);
                catch
                    lat(jj)     = NaN;
                    lon(jj)     = NaN;
                    amp(jj)     = NaN;
                    phs(jj)     = NaN;
                end
            end
            if xx == 1; hscatter2   = scatter(lon, lat, 20, amp, 'filled'); end
            if xx == 2; hscatter2   = scatter(lon, lat, 20, phs, 'filled'); end
            set(hscatter, 'MarkerEdgeColor', 'w');
            set(hscatter2, 'MarkerEdgeColor', 'm');
            axis equal
            xlim([-205 -135]);
            ylim([48 81]);
            if xx == 1; caxis([0 numbers.high95]); end
            if xx == 2; caxis([0 360]); cmap = cmocean('phase'); colormap(cmap); end
            xlabel('longitude [\circ]');
            ylabel('latitude [\circ]');
            if xx== 1; hc = colorbar; ylabel(hc, [components_wanted{ii}, ' amplitude [m]']); end
            if xx== 2; hc = colorbar; ylabel(hc, [components_wanted{ii}, ' phase [degrees]']); end
            if xx==1; figure_png = [components_wanted{ii}, '_amplitude.png']; end
            if xx==2; figure_png = [components_wanted{ii}, '_phase.png']; end
            legend([h hscatter(1) hscatter2], 'D-Flow FM', 'FES2014', 'observed', 'location', 'northwest');
            print('-dpng','-r300',figure_png); close all
        end
    end
end