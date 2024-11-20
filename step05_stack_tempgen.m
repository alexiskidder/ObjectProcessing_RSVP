%% stack tmep gen results and average
% Nove 22 2023 - AK 


%% set up

close all
clear all

    %addpath('~/CoSMoMVPA/mvpa/') % cosmomvpa

rootdir = "/Users/f004p7b/Documents/australia/objects_project/analysis_bids/data/derivatives/tempgen/";
outdir = "/Users/f004p7b/Documents/australia/objects_project/plots/tempgen/" ;
parts = dir("/Users/f004p7b/Documents/australia/objects_project/analysis_bids/data/derivatives/tempgen/*anim_tempgen_orig.mat");
dim_values = load(rootdir + "dim_values.mat");
dim_labels = load(rootdir + "dim_labels.mat");

% where to put the files 
arorig_all = [];
arbin_all = [];

animorig_all = [];
animbin_all = [];

catorig_all = [];
catbin_all = []; 

imorig_all = [];
imbin_all = [];

%% get appropriate files for each participant 

for g = 1:length(parts)

    basefil = parts(g).name;
    subj = basefil(1:6);
    disp("stacking " + subj)

    arfile = load(rootdir + subj + "_asprat_all.mat");
    animfile = load(rootdir + subj + "tempgen_anim_all.mat");
    catfile = load(rootdir + subj + "tempgen_cat_all.mat"); 
    imfile = load(rootdir + subj + "tempgen_image_all.mat");

    for mid = 1:4

        if mid == 1 
            file = arfile;  
            origds = file.allaspr_out{1};
            binds  = file.allaspr_out{2};

        elseif mid == 2 
            file = animfile;
            origds = file.allanim_out{1};
            binds  = file.allanim_out{2};

        elseif mid == 3
            file = catfile; 
            origds = file.allcat_out{1};
            binds  = file.allcat_out{2};

        else
            file = imfile; 
            origds = file.allimage_out{1};
            binds  = file.allimage_out{2};

        end

    ds_o = cosmo_stack(origds);
    ds_o = cosmo_dim_transpose(ds_o,'train_time');
    ds_o = cosmo_dim_transpose(ds_o,'test_time');
    ds_o = cosmo_fx(ds_o,@mean);
    [Xo,dim_labels,dim_values] = cosmo_unflatten(ds_o);

    %origtempout{g,1} = Xo;

    ds_b = cosmo_stack(binds);
    ds_b = cosmo_dim_transpose(ds_b,'train_time');
    ds_b = cosmo_dim_transpose(ds_b,'test_time');
    ds_b = cosmo_fx(ds_b,@mean);
    [Xb,dim_labels,dim_values] = cosmo_unflatten(ds_b);

    %bintempout{g,1} = Xb;

    if mid == 1 
        arorig_all{g,1} = Xo;
        arbin_all{g,1} = Xb;
    elseif mid == 2
        animorig_all{g,1} = Xo;
        animbin_all{g,1} = Xb;
    elseif mid == 3
        catorig_all{g,1} = Xo;
        catbin_all{g,1} = Xb;
    else
        imorig_all{g,1} = Xo;
        imbin_all{g,1} = Xb;

    end

end

end


%% save files 

orig_ar_outfile = rootdir + "all_p=" + length(parts) + "_orig_asprat.mat";
save(orig_ar_outfile,"arorig_all");

bin_ar_outfile = rootdir + "all_p=" + length(parts) + "_bin_asprat.mat";
save(bin_ar_outfile,"arbin_all");

orig_anim_outfile = rootdir + "all_p=" + length(parts) + "_orig_anim.mat";
save(orig_anim_outfile,"animorig_all");

bin_anim_outfile = rootdir + "all_p=" + length(parts) + "_bin_anim.mat";
save(bin_anim_outfile,"animbin_all");

orig_cat_outfile = rootdir + "all_p=" + length(parts) + "_orig_cat.mat";
save(orig_cat_outfile,"catorig_all");

bin_cat_outfile = rootdir + "all_p=" + length(parts) + "_bin_cat.mat";
save(bin_cat_outfile,"catbin_all");

orig_im_outfile = rootdir + "all_p=" + length(parts) + "_orig_im.mat";
save(orig_im_outfile,"imorig_all");

bin_im_outfile = rootdir + "all_p=" + length(parts) + "_bin_im.mat";
save(bin_im_outfile,"imbin_all");


%% average together 
% coded the above portion slightly wrong, so can just use the last cell in
% each combined struct

%%% original - aspect ratio
    arorigmat = cell2mat(arorig_all);
    arorigmat = permute(arorigmat, [2 3 1]);% want participants to be the 3rd dimension 
    mean_arorigall = mean(arorigmat,3);

% binarized - aspect ratio
    arbinmat = cell2mat(arbin_all);
    arbinmat = permute(arbinmat, [2 3 1]);% want participants to be the 3rd dimension 
    mean_arbinall = mean(arbinmat,3);

%%% original - animacy
    animorigmat = cell2mat(animorig_all);
    animorigmat = permute(animorigmat, [2 3 1]);% want participants to be the 3rd dimension 
    mean_animorigall = mean(animorigmat,3);

% binarized - animacy
    animbinmat = cell2mat(animbin_all);
    animbinmat = permute(animbinmat, [2 3 1]);% want participants to be the 3rd dimension 
    mean_animbinall = mean(animbinmat,3);

%%% original - category
    catorigmat = cell2mat(catorig_all);
    catorigmat = permute(catorigmat, [2 3 1]);% want participants to be the 3rd dimension 
    mean_catorigall = mean(catorigmat,3);

% binarized - category
    catbinmat = cell2mat(catbin_all);
    catbinmat = permute(catbinmat, [2 3 1]);% want participants to be the 3rd dimension 
    mean_catbinall = mean(catbinmat,3);

%%% original - image
    imorigmat = cell2mat(imorig_all);
    imorigmat = permute(imorigmat, [2 3 1]);% want participants to be the 3rd dimension 
    mean_imorigall = mean(imorigmat,3);

% binarized - image
    imbinmat = cell2mat(imbin_all);
    imbinmat = permute(imbinmat, [2 3 1]);% want participants to be the 3rd dimension 
    mean_imbinall = mean(imbinmat,3);

%% find the means for each the color maps!

ar_orig_low_high = [];
anim_orig_low_high = [];
cat_orig_low_high = [];
im_orig_low_high = [];

temp_max = max(mean_arorigall);
temp_min = min(mean_arorigall);
current_max = max(temp_max);
current_min = min(temp_min);
ar_orig_low_high(1,1) = current_min;
ar_orig_low_high(1,2) = current_max;

temp_max = max(mean_arbinall);
temp_min = min(mean_arbinall);
current_max = max(temp_max);
current_min = min(temp_min);
ar_orig_low_high(2,1) = current_min;
ar_orig_low_high(2,2) = current_max;

temp_max = max(mean_animorigall);
temp_min = min(mean_animorigall);
current_max = max(temp_max);
current_min = min(temp_min);
anim_orig_low_high(1,1) = current_min;
anim_orig_low_high(1,2) = current_max;

temp_max = max(mean_animbinall);
temp_min = min(mean_animbinall);
current_max = max(temp_max);
current_min = min(temp_min);
anim_orig_low_high(2,1) = current_min;
anim_orig_low_high(2,2) = current_max;

temp_max = max(mean_catorigall);
temp_min = min(mean_catorigall);
current_max = max(temp_max);
current_min = min(temp_min);
cat_orig_low_high(1,1) = current_min;
cat_orig_low_high(1,2) = current_max;

temp_max = max(mean_catbinall);
temp_min = min(mean_catbinall);
current_max = max(temp_max);
current_min = min(temp_min);
cat_orig_low_high(2,1) = current_min;
cat_orig_low_high(2,2) = current_max;

temp_max = max(mean_imorigall);
temp_min = min(mean_imorigall);
current_max = max(temp_max);
current_min = min(temp_min);
im_orig_low_high(1,1) = current_min;
im_orig_low_high(1,2) = current_max;

temp_max = max(mean_imbinall);
temp_min = min(mean_imbinall);
current_max = max(temp_max);
current_min = min(temp_min);
im_orig_low_high(2,1) = current_min;
im_orig_low_high(2,2) = current_max;


%% plot 

%% custom color maps 

% create color map - teal to orange

lines_map = lines;

% Define the start, middle, and end colors
start_color = lines_map(1, :);  % First color in 'lines' colormap (teal)
end_color = lines_map(2, :);    % Second color in 'lines' colormap (orange)
white_color = [1 1 1];          % White color

% Number of intermediate colors
num_steps = 100;

% Generate a smooth gradient from start color to white to end color
step1 = linspace(start_color(1), white_color(1), num_steps/2);
step2 = linspace(white_color(1), end_color(1), num_steps/2);
comb1 = [step1,step2];

step3 = linspace(start_color(2), white_color(2), num_steps/2);
step4 = linspace(white_color(2), end_color(2), num_steps/2);
comb2 = [step3,step4];

start_to_white = [linspace(start_color(1), white_color(1), num_steps/2), ...
                  linspace(white_color(1), end_color(1), num_steps/2)];
              
start_to_white = [linspace(start_color(2), white_color(2), num_steps/2), ...
                  linspace(white_color(2), end_color(2), num_steps/2)];

start_to_white = [linspace(start_color(3), white_color(3), num_steps/2), ...
                  linspace(white_color(3), end_color(3), num_steps/2)];

% Combine the RGB components to create the colormap
custom_map = [comb1', comb2', start_to_white'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Define the start and end colors
start_color = [0.6 0.8 1];    % Light blue color
end_color = [0.2, 0, 0.4];
%end_color = parula(1);        % First color in 'parula' colormap

% Number of intermediate colors
num_steps = 100;

% Generate a smooth gradient from start color to end color
r = linspace(start_color(1), end_color(1), num_steps);
g = linspace(start_color(2), end_color(2), num_steps);
b = linspace(start_color(3), end_color(3), num_steps);

% Combine the RGB components to create the colormap
blues = [r', g', b'];

% Plot a colorbar to visualize the colormap
imagesc(1:num_steps, [0 1], reshape(blues, [num_steps 1 3]));
colormap(blues);
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Define the start, middle, and end colors
start_color = [0.2, 0.2, 0.2];  % Dark grey color
%end_color = [0, 0.2, 0.4];      % Desired end color
end_color = [0, 0.3, 0.4];      % Desired end color
white_color = [1, 1, 1];        % White color

% Number of intermediate colors
num_steps = 100;

% Generate a smooth gradient from start color to white to end color
stepa = linspace(start_color(1), white_color(1), num_steps/2);
stepb = linspace(white_color(1), end_color(1), num_steps/2);
comba = [stepa, stepb];

stepc = linspace(start_color(2), white_color(2), num_steps/2);
stepd = linspace(white_color(2), end_color(2), num_steps/2);
combb = [stepc, stepd];

start_to_white = [linspace(start_color(3), white_color(3), num_steps/2), ...
                  linspace(white_color(3), end_color(3), num_steps/2)];

% Combine the RGB components to create the colormap
grey_blue_map = [comba', combb', start_to_white'];

% Plot a colorbar to visualize the colormap
imagesc(1:num_steps, [0 1], reshape(grey_blue_map, [num_steps 1 3]));
colormap(grey_blue_map);
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Define the start and end colors
start_color = [0.7529, 0.8471, 0.9176];    % Light white blue color
end_color = [0, 0.3, 0.4];  %teal like color

% Number of intermediate colors
num_steps = 100;

% Generate a smooth gradient from start color to end color
r = linspace(start_color(1), end_color(1), num_steps);
g = linspace(start_color(2), end_color(2), num_steps);
b = linspace(start_color(3), end_color(3), num_steps);

% Combine the RGB components to create the colormap
new_blues = [r', g', b'];

% Plot a colorbar to visualize the colormap
imagesc(1:num_steps, [0 1], reshape(new_blues, [num_steps 1 3]));
colormap(new_blues);
colorbar;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Define the start and end colors
start_color = [0.7490, 0.6078, 0.2275]; % mustard yellow
end_color = [0, 0.3, 0.4];  %teal like color

% Number of intermediate colors
num_steps = 100;

% Generate a smooth gradient from start color to end color
r = linspace(start_color(1), end_color(1), num_steps);
g = linspace(start_color(2), end_color(2), num_steps);
b = linspace(start_color(3), end_color(3), num_steps);

% Combine the RGB components to create the colormap
yel_blue = [r', g', b'];

% Plot a colorbar to visualize the colormap
imagesc(1:num_steps, [0 1], reshape(yel_blue, [num_steps 1 3]));
colormap(yel_blue);
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the colors
col1 = [0.4, 0, 0.6]; %purple
col2 = [1, 0.8, 0.8];%lightpink
col3 = [0.5, 0.5, 0.9]; %blue
col4 = [0, 0.5, 0.5];%teal

% Number of steps
num_steps = 300;

% Generate the colormap
r = [linspace(col1(1), col2(1), num_steps/3), ...
     linspace(col2(1), col3(1), num_steps/3), ...
     linspace(col3(1), col4(1), num_steps/3)];
g = [linspace(col1(2), col2(2), num_steps/3), ...
     linspace(col2(2), col3(2), num_steps/3), ...
     linspace(col3(2), col4(2), num_steps/3)];
b = [linspace(col1(3), col2(3), num_steps/3), ...
     linspace(col2(3), col3(3), num_steps/3), ...
     linspace(col3(3), col4(3), num_steps/3)];

% Combine the RGB components to create the colormap
yelmulti_map = [r', g', b'];

% Plot a colorbar to visualize the colormap
imagesc(1:num_steps, [0 1], reshape(yelmulti_map, [num_steps 1 3]));
colormap(yelmulti_map);
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the colors
col1 = [0.5, 0, 0.5]; %purple
col2 = [1, 0.8, 0.8];%lightpink
col3 = [0, 0, 1]; %blue
col4 = [0, 0.5, 0];%green 

% Number of steps
num_steps = 300;

% Generate the colormap
r = [linspace(col1(1), col2(1), num_steps/3), ...
     linspace(col2(1), col3(1), num_steps/3), ...
     linspace(col3(1), col4(1), num_steps/3)];
g = [linspace(col1(2), col2(2), num_steps/3), ...
     linspace(col2(2), col3(2), num_steps/3), ...
     linspace(col3(2), col4(2), num_steps/3)];
b = [linspace(col1(3), col2(3), num_steps/3), ...
     linspace(col2(3), col3(3), num_steps/3), ...
     linspace(col3(3), col4(3), num_steps/3)];

% Combine the RGB components to create the colormap
multi_map = [r', g', b'];

% Plot a colorbar to visualize the colormap
imagesc(1:num_steps, [0 1], reshape(multi_map, [num_steps 1 3]));
colormap(multi_map);
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the colors
col1 = [0, 0.5, 0.5];%teal
col2 = [1, 0.8, 0.8];%lightpink
col3 = [0.5, 0.5, 0.9]; %blue
col4 = [0.4, 0, 0.6]; %purple

% Number of steps
num_steps = 300;

% Generate the colormap
r = [linspace(col1(1), col2(1), num_steps/3), ...
     linspace(col2(1), col3(1), num_steps/3), ...
     linspace(col3(1), col4(1), num_steps/3)];
g = [linspace(col1(2), col2(2), num_steps/3), ...
     linspace(col2(2), col3(2), num_steps/3), ...
     linspace(col3(2), col4(2), num_steps/3)];
b = [linspace(col1(3), col2(3), num_steps/3), ...
     linspace(col2(3), col3(3), num_steps/3), ...
     linspace(col3(3), col4(3), num_steps/3)];

% Combine the RGB components to create the colormap
tealmulti_map = [r', g', b'];

% Plot a colorbar to visualize the colormap
imagesc(1:num_steps, [0 1], reshape(tealmulti_map, [num_steps 1 3]));
colormap(tealmulti_map);
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the colors
col1 = [0, 0, 0.4];%blue
col2 = [1, 0.8, 0.8];%lightpink
col3 = [0.8, 0.4, 0]; %orange
col4 = [0.769, 0.686, 0.286]; %yellow

% Number of steps
num_steps = 300;

% Generate the colormap
r = [linspace(col1(1), col2(1), num_steps/3), ...
     linspace(col2(1), col3(1), num_steps/3), ...
     linspace(col3(1), col4(1), num_steps/3)];
g = [linspace(col1(2), col2(2), num_steps/3), ...
     linspace(col2(2), col3(2), num_steps/3), ...
     linspace(col3(2), col4(2), num_steps/3)];
b = [linspace(col1(3), col2(3), num_steps/3), ...
     linspace(col2(3), col3(3), num_steps/3), ...
     linspace(col3(3), col4(3), num_steps/3)];

% Combine the RGB components to create the colormap
sun_map = [r', g', b'];

% Plot a colorbar to visualize the colormap
imagesc(1:num_steps, [0 1], reshape(sun_map, [num_steps 1 3]));
colormap(sun_map);
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the colors
col1 = [0, 0, .8];%blue
col2 = [1, 0.8, 0.8];%lightpink
col3 = [0.8, 0.4, 0]; %orange
col4 = [0.769, 0.686, 0.286]; %yellow

% Number of steps
num_steps = 300;

% Generate the colormap
r = [linspace(col1(1), col2(1), num_steps/3), ...
     linspace(col2(1), col3(1), num_steps/3), ...
     linspace(col3(1), col4(1), num_steps/3)];
g = [linspace(col1(2), col2(2), num_steps/3), ...
     linspace(col2(2), col3(2), num_steps/3), ...
     linspace(col3(2), col4(2), num_steps/3)];
b = [linspace(col1(3), col2(3), num_steps/3), ...
     linspace(col2(3), col3(3), num_steps/3), ...
     linspace(col3(3), col4(3), num_steps/3)];

% Combine the RGB components to create the colormap
sun2_map = [r', g', b'];

% Plot a colorbar to visualize the colormap
imagesc(1:num_steps, [0 1], reshape(sun2_map, [num_steps 1 3]));
colormap(sun2_map);
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the colors
col1 = [0, 0, .8];%blue
col2 = [1, 0.8, 0.8];%lightpink
col3 = [0.8, 0.4, 0]; %orange
col4 = [0.5, 0, 0]; %maroon

% Number of steps
num_steps = 300;

% Generate the colormap
r = [linspace(col1(1), col2(1), num_steps/3), ...
     linspace(col2(1), col3(1), num_steps/3), ...
     linspace(col3(1), col4(1), num_steps/3)];
g = [linspace(col1(2), col2(2), num_steps/3), ...
     linspace(col2(2), col3(2), num_steps/3), ...
     linspace(col3(2), col4(2), num_steps/3)];
b = [linspace(col1(3), col2(3), num_steps/3), ...
     linspace(col2(3), col3(3), num_steps/3), ...
     linspace(col3(3), col4(3), num_steps/3)];

% Combine the RGB components to create the colormap
sun3_map = [r', g', b'];

% Plot a colorbar to visualize the colormap
imagesc(1:num_steps, [0 1], reshape(sun3_map, [num_steps 1 3]));
colormap(sun3_map);
colorbar;



%% back to plotting 

%colormap(custom_map)
%colormap(blues)
%colormap(grey_blue_map)
%colormap(new_blues)

% aspect ratio 
% original
    colormap(sun3_map)
    figure(1);clf
    imagesc(dim_values{1},dim_values{2},squeeze(mean_arorigall));hold on
    axis square
    xlabel(strrep(dim_labels{1},'_',' '))
    ylabel(strrep(dim_labels{2},'_',' '))
    a=gca;
    a.YDir='normal';
    plot(a.XLim,a.YLim,'--','Color',[0.35 0.35 0.35])
    %c = colorbar();c.Label.String='accuracy'
    %c = colorbar(); caxis([min(ar_orig_low_high(:,1)), max(ar_orig_low_high(:,2))])
    c = colorbar(); caxis([.44, .63])
    %title('aspect ratio - original stimuli');
    filename = outdir + "all_p=" + length(parts) + "_aspect_orig_sun3.png";
    set(gca,'FontSize',16)
    set(gcf,'color','w')
    Image = getframe(gcf) ; % gets info from the opened figure
    imwrite(Image.cdata, filename);

% binarized
    colormap(sun3_map)
    figure(1);clf
    imagesc(dim_values{1},dim_values{2},squeeze(mean_arbinall));hold on
    axis square
    xlabel(strrep(dim_labels{1},'_',' '))
    ylabel(strrep(dim_labels{2},'_',' '))
    a=gca;
    a.YDir='normal';
    plot(a.XLim,a.YLim,'--','Color',[0.35 0.35 0.35])
    %c = colorbar();c.Label.String='accuracy'
    %c = colorbar(); caxis([min(ar_orig_low_high(:,1)), max(ar_orig_low_high(:,2))])
    c = colorbar(); caxis([.44, .63])
    %title('aspect ratio - binarized stimuli');
    filename = outdir + "all_p=" + length(parts) + "_aspect_bin_sun3.png";
    set(gca,'FontSize',16)
    set(gcf,'color','w')
    Image = getframe(gcf) ; % gets info from the opened figure
    imwrite(Image.cdata, filename);

%%% animacy

% original
    colormap(sun3_map)
    figure(1);clf
    imagesc(dim_values{1},dim_values{2},squeeze(mean_animorigall));hold on
    axis square
    xlabel(strrep(dim_labels{1},'_',' '))
    ylabel(strrep(dim_labels{2},'_',' '))
    a=gca;
    a.YDir='normal';
    plot(a.XLim,a.YLim,'--','Color',[0.35 0.35 0.35])
    %c = colorbar();c.Label.String='accuracy'
    %c = colorbar(); caxis([min(anim_orig_low_high(:,1)), max(anim_orig_low_high(:,2))])
    c = colorbar(); caxis([.44, .63])
    %title('animacy - original stimuli');
    filename = outdir + "all_p=" + length(parts) + "_animacy_orig_sun3.png";
    set(gca,'FontSize',16)
    set(gcf,'color','w')
    Image = getframe(gcf) ; % gets info from the opened figure
    imwrite(Image.cdata, filename);

% binarized
    colormap(sun3_map)
    figure(1);clf
    imagesc(dim_values{1},dim_values{2},squeeze(mean_animbinall));hold on
    axis square
    xlabel(strrep(dim_labels{1},'_',' '))
    ylabel(strrep(dim_labels{2},'_',' '))
    a=gca;
    a.YDir='normal';
    plot(a.XLim,a.YLim,'--','Color',[0.35 0.35 0.35])
    c = colorbar(); caxis([.44, .63])
    %c = colorbar();c.Label.String='accuracy'
    %c = colorbar(); caxis([min(anim_orig_low_high(:,1)), max(anim_orig_low_high(:,2))])
    %title('animacy - binarized stimuli');
    filename = outdir + "all_p=" + length(parts) + "_animacy_bin_sun3.png";
    set(gca,'FontSize',16)
    set(gcf,'color','w')
    Image = getframe(gcf) ; % gets info from the opened figure
    imwrite(Image.cdata, filename);


%%% category

% original
    colormap(sun3_map)
    figure(1);clf
    imagesc(dim_values{1},dim_values{2},squeeze(mean_catorigall));hold on
    axis square
    xlabel(strrep(dim_labels{1},'_',' '))
    ylabel(strrep(dim_labels{2},'_',' '))
    a=gca;
    a.YDir='normal';
    plot(a.XLim,a.YLim,'--','Color',[0.35 0.35 0.35])
    %c = colorbar();c.Label.String='accuracy'
    c = colorbar(); caxis([min(cat_orig_low_high(:,1)), max(cat_orig_low_high(:,2))])
    %title('category - original stimuli');
    filename = outdir + "all_p=" + length(parts) + "_category_orig_sun3.png";
    set(gca,'FontSize',16)
    set(gcf,'color','w')
    Image = getframe(gcf) ; % gets info from the opened figure
    imwrite(Image.cdata, filename);

% binarized
    colormap(sun3_map)
    figure(1);clf
    imagesc(dim_values{1},dim_values{2},squeeze(mean_catbinall));hold on
    axis square
    xlabel(strrep(dim_labels{1},'_',' '))
    ylabel(strrep(dim_labels{2},'_',' '))
    a=gca;
    a.YDir='normal';
    plot(a.XLim,a.YLim,'--','Color',[0.35 0.35 0.35])
    %c = colorbar();c.Label.String='accuracy'
    c = colorbar(); caxis([min(cat_orig_low_high(:,1)), max(cat_orig_low_high(:,2))])
    %title('category - binarized stimuli');
    filename = outdir + "all_p=" + length(parts) + "_category_bin_sun3.png";
    set(gca,'FontSize',16)
    set(gcf,'color','w')
    Image = getframe(gcf) ; % gets info from the opened figure
    imwrite(Image.cdata, filename);























