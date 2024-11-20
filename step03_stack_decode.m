%% shortened decode stack
% November 8, 2023 - AK 

%% set-up & load

close all;
clear all;
rootdir = '/Users/f004p7b/Documents/australia/objects_project/analysis_bids/data/' ;
decodir = rootdir + "derivatives/decoding/";
cd(decodir) 
files = dir('*res_decoding.mat'); %get the files that correspond to all participants in the folder

%cosmo and fieldtrip info (don't know if I need this or not)
    addpath('~/CoSMoMVPA/mvpa/')
    addpath('~/fieldtrip-20231015');
    ft_defaults;


%% stack the participants decoding results 

% location for information from all participants - original
animacy_all = [];
cat_all = [];
asp_all = [];
image_all = [];

% location for information from all participants - binarized

bin_animacy_all = [];
bin_cat_all = [];
bin_asp_all = [];
bin_image_all = [];

% make combined file/info

for pnum = 1:length(files)

    fnam = files(pnum).name;
    pfile = decodir + fnam;
    
    k = load(pfile);
    
    animacy_all{pnum,1} = k.res_cell_anim{1, 1};
    cat_all{pnum,1} = k.res_cell_cat{1,1};
    asp_all{pnum,1} = k.res_cell_asp{1,1};
    image_all{pnum,1} = k.res_cell_image{1,1};

     bin_animacy_all{pnum,1} = k.res_cell_anim{1, 2};
     bin_cat_all{pnum,1} = k.res_cell_cat{1,2};
     bin_asp_all{pnum,1} = k.res_cell_asp{1,2};
     bin_image_all{pnum,1} = k.res_cell_image{1,2};
end

% stack

res_animacy = cosmo_stack(animacy_all);
res_cat = cosmo_stack(cat_all);
res_asp = cosmo_stack(asp_all);
res_image = cosmo_stack(image_all);

res_bin_animacy = cosmo_stack(bin_animacy_all);
res_bin_cat = cosmo_stack(bin_cat_all);
res_bin_asp = cosmo_stack(bin_asp_all);
res_bin_image = cosmo_stack(bin_image_all);

%% save file

% get the number of participants for the dataset

lastpar = files((length(files)),1).name;
dsnum = lastpar(5:6);

% save file

filnam = decodir + "short_origstim-decoding-allthrough" + dsnum + ".mat";
binfilnam = decodir + "short_binstim-decoding-allthrough" + dsnum + ".mat";

save(filnam,"res_image","res_animacy","res_cat", "res_asp");
save(binfilnam,"res_bin_animacy", "res_bin_cat", "res_bin_asp", "res_bin_image");

%% average

animacy_av{1} = mean(res_animacy.samples,1);
animacy_av{2} = mean(res_bin_animacy.samples,1);

aspect_av{1} = mean(res_asp.samples,1);
aspect_av{2} = mean(res_bin_asp.samples,1);

category_av{1} = mean(res_cat.samples,1);
category_av{2} = mean(res_bin_cat.samples,1);

image_av{1} = mean(res_image.samples,1);
image_av{2} = mean(res_bin_image.samples,1);


%% save file

annam =  decodir + "short_avg-animacydecoding-allthrough" + dsnum + ".mat";
save(annam,"animacy_av");

aspnam =  decodir + "short_avg-aspectratiodecoding-allthrough" + dsnum + ".mat";
save(aspnam,"aspect_av");

catnam =  decodir + "short_avg-categorydecoding-allthrough" + dsnum + ".mat";
save(catnam,"category_av");

imanam = decodir + "short_avg-imagedecoding-allthrough" + dsnum + ".mat";
save(imanam,"image_av");


%% plot
figure(1);clf

% image
subplot(4,1,1)
timelab = res_image.a.fdim.values{1};
for p = 1:2
    plot(timelab,image_av{p}); hold on
end
plot(timelab,1/52+0*timelab,'k--')
legend({'original','binarized','chance'})
xlim(minmax(timelab))
title("decoding image - average of " + length(files) + " participants")

% category
subplot(4,1,2)
timelab = res_cat.a.fdim.values{1};
for p = 1:2
    plot(timelab,category_av{p}); hold on
end
plot(timelab,1/4+0*timelab,'k--')
legend({'original','binarized','chance'})
xlim(minmax(timelab))
title("decoding category - average of " + length(files) + " participants")

% aspect ratio
subplot(4,1,3)
timelab = res_asp.a.fdim.values{1};
for p = 1:2
    plot(timelab,aspect_av{p}); hold on
end
plot(timelab,1/2+0*timelab,'k--')
legend({'original','binarized','chance'})
xlim(minmax(timelab))
title("decoding aspect ratio - average of " + length(files) + " participants")

% animacy
subplot(4,1,4)
timelab = res_animacy.a.fdim.values{1};
for p = 1:2
    plot(timelab,animacy_av{p}); hold on
end

plot(timelab,1/2+0*timelab,'k--')
legend({'original','binarized','chance'})
xlim(minmax(timelab))
title("decoding animacy - average of " + length(files) + " participants")


%% save plot

loc = '/Users/f004p7b/Documents/australia/objects_project/plots/decoding_figures/';
tempname = "short_avg_n" + length(files) + "_decoding.png";
tn = loc + tempname;
print(gcf,'-dpng','-r500',tn)
im=imread([tn]);
[i,j]=find(mean(im,3)<255);margin=0;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[tn '.png'],'png');
 % note: imwrite sometimes throws a weird error, but it still saves and
 % everything is fine so ignore it 













