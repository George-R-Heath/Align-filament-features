
min_ccr =0.2;      %minimum correlation coefficient value 
use_ref = 0;       %Reference frame selection 
%Set to 0 go to line 15 to type in an image region
%Set to 1 uses correlation average (only when re-running code)
%Set to 2 to load new image as reference
%Set to 3 to keep current reference (only when re-running code)

sort_stack = 0;    %set to 1 to sort the stack by correlation coefficient
fr = 0.5;          %if sorting stack set cut off correlation coefficient value

%%
stitch = [];
for i = 1:numel(Fco_saved)
stitch = [stitch, dig_straight{i}]; 
end

%set use_ref = 3 to keep current ref (only if repeating runs)
if use_ref ==0  
ref = stitch(:,51:100); %if use_ref = 0 enter pixel region of filament as reference

elseif use_ref ==1
ref = avg_clip;         %if use_ref = 1 average cc will be used as reference (when running the code to refine)

elseif use_ref ==2              
[f,path] = uigetfile('*.tif');   %if use_ref = 2 to open a new file as a reference 
f = fullfile(path,f);
info = imfinfo(f);
d1 = imread(f,1);
d1 = im2double(d1);
ref = d1;
end
paddingy = 20;%pad(1)*2;
paddingx = 20;%pad(2)*2;

ref_sd = size(ref);
clear I cent_co cent_co_top srr sort_ssr clip ccr x y

pad = size(ref);
min_peak_sep = pad(2)/3;

if (paddingy~=0)
    stitch(end+paddingy,:)=0;
    stitch(:,end+paddingx)=0;
   stitch = imtranslate(stitch,[paddingx/2,paddingy/2]);
end

edg = 3;

ccr = normxcorr2(ref,stitch);
ccr(1:paddingy/2,:) = 0;
ccr(end-paddingy/2:end,:) = 0;
ccr(:,1:paddingx/2) = 0;
ccr(:,end-paddingx/2:end) = 0;

ccr = imgaussfilt(ccr,1);

 sd=size(ccr);
                [x y]=find(ccr(edg:sd(1)-edg,edg:sd(2)-edg));
                
                % initialize outputs
                cent=[]; cent_h=[];
                cent_map=zeros(sd);
                
                x=x+edg-1;
                y=y+edg-1;
                for j=1:length(y)
                    if (ccr(x(j),y(j))>=ccr(x(j)-1,y(j)-1 )) &&...
                            (ccr(x(j),y(j))>ccr(x(j)-1,y(j))) &&...
                            (ccr(x(j),y(j))>=ccr(x(j)-1,y(j)+1)) &&...
                            (ccr(x(j),y(j))>ccr(x(j),y(j)-1)) && ...
                            (ccr(x(j),y(j))>ccr(x(j),y(j)+1)) && ...
                            (ccr(x(j),y(j))>=ccr(x(j)+1,y(j)-1)) && ...
                            (ccr(x(j),y(j))>ccr(x(j)+1,y(j))) && ... 
                              (ccr(x(j),y(j))>min_ccr) && ...
                            (ccr(x(j),y(j))>=ccr(x(j)+1,y(j)+1));

                        cent = [cent ;  y(j) ; x(j)];
                        cent_map(x(j),y(j))=cent_map(x(j),y(j))+1; % if a binary matrix output is desired
                    end
                end
cent_co = [cent(1:2:end),cent(2:2:end)];
for i=1:numel(cent_co(:,1))
cent_h(i) = ccr(cent_co(i,2),cent_co(i,1));
end

 keep = [];
    %remove overlapping correlation neighbours
    if any(cent_co)
    for j = 1:numel(cent_co(:,1))
    [ne,d] = knnsearch(cent_co,cent_co(j,:),'k',5);
     pos = d<min_peak_sep;
        if sum(pos)>1
        [del_pos del_h] = max(cent_h(ne(pos)));
        keep(j) = (ne(del_h));
        else
        keep(j) = j;
        end
    end
keep = unique(keep'.');
cent_co = cent_co(keep,:);
cent_h = cent_h(:,keep);
     end

[sort_ssr, I] = sort(cent_h,'descend');
pos = sort_ssr>fr; top = sum(pos);
cent_co_top = cent_co(I(1:top),:);

figure(2)
plot(sort_ssr)
hold on
plot(sort_ssr(1:top),'o')

figure(4)
tiledlayout(3,1, 'Padding', 'none', 'TileSpacing', 'none');  
nexttile
imagesc(ccr)
title('Cross Correlation')
hold on
plot(cent_co(:,1),cent_co(:,2),'.')
plot(cent_co_top(:,1),cent_co_top(:,2),'o')
colormap(lutafm)
nexttile
plot(cent_co(:,1),cent_h,'-o')
grid on
xlim([0 sd(2)])
nexttile
imagesc(stitch)
xlim([-ref_sd(2)/2 sd(2)])
set(gca,'YDir','normal')

figure('Position',[10 100 1000 200*ceil(numel(cent_h)/10)])
if sort_stack == 1   
tiledlayout(ceil((top+2)/10),10, 'Padding', 'none', 'TileSpacing', 'none');
else
    tiledlayout(ceil((numel(cent_h)+2)/10),10, 'Padding', 'none', 'TileSpacing', 'none');
end
nexttile
    imagesc(ref)
    title('Reference')
    set(gca,'YDir','normal')
    colormap(lutafm)
      avg_clip = zeros(ref_sd(1), ref_sd(2));
if sort_stack == 1
    for i=1:top
    x = round(cent_co_top(i,1)-ref_sd(2)/2);
    y = round(cent_co_top(i,2)-ref_sd(1)/2);
        clip{i} = stitch((y-round(ref_sd(1)/2)+1):(y+round(ref_sd(1)/2)), (x-round(ref_sd(2)/2)+1):x+round(ref_sd(2)/2));
        avg_clip = avg_clip + clip{i};
nexttile
    imagesc(clip{i})
    title(i)
    set(gca,'YDir','normal')
    end
    avg_clip = avg_clip/numel(top);
    nexttile
    imagesc(avg_clip)
    title('CC AVG')
    set(gca,'YDir','normal')
else
    for i=1:numel(cent_h)
    x = round(cent_co(i,1)-ref_sd(2)/2);
    y = round(cent_co(i,2)-ref_sd(1)/2);
    clip{i} = stitch((y-round(ref_sd(1)/2)+1):(y+round(ref_sd(1)/2)), (x-round(ref_sd(2)/2)+1):x+round(ref_sd(2)/2));
    avg_clip = avg_clip + clip{i};
    nexttile
    imagesc(clip{i})
    title(i)
    set(gca,'YDir','normal')
    end  
    avg_clip = avg_clip/numel(cent_h);
    nexttile
    imagesc(avg_clip)
    title('CC AVG')
    set(gca,'YDir','normal')
end


data = uint32(flip(clip{1}*1000,1));
outputFileName = 'cc stack.tif';

t = Tiff(outputFileName,'w');
tagstruct.ImageLength     = size(data,1);
tagstruct.ImageWidth      = size(data,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 32;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct)
t.write(data);
t.close();

for i = 2:numel(clip)
data = uint32(flip(clip{i}*1000,1));
t = Tiff(outputFileName,'a');
tagstruct.ImageLength     = size(data,1);
tagstruct.ImageWidth      = size(data,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 32;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct)
t.write(data);
t.close();
end
    