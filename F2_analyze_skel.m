%pixel interpolation, digital straighten

prof_smoo = 20;   %set degree of filament profile smoothing 
pleng = 50;       %set cross section length
show_cross=0;     %set to 1 to see direction of cross sections

minlength = 15;
%%
clear perp
loctor = 0.2; % set the location tolerance (1=no tolerance) away from centre of peak postion for width measure
F = griddedInterpolant(double(A)); %to interpolate pixel values in digital straighten 
F.Method = 'cubic';
wd_analysis = 1;  

          for i=1:nclus
            [r c]=find(clusters==i);
             if length(r)<minlength
                
            else
               
 clustersedit=clusters;
                idx=clustersedit~=i;
                clustersedit(idx) = 0 ;
                idx=clustersedit==i;
               % maskclus(c1:c2,c1:c2)=dataout.*(idx);
                %raw_data2_clus=single(maskclus(c1:c2,c1:c2));
                
                %% Skeletonising
                clustersedit = bwmorph(clustersedit,'bridge');
                clustersedit = bwmorph(clustersedit,'clean');
                clustersedit = bwmorph(clustersedit,'close');
                clustersedit = bwmorph(clustersedit,'fill');
                clustersedit=bwmorph(clustersedit, 'thin', Inf) ;
                clustersedit = bwskel(logical(clustersedit),'MinBranchLength',mbl);
                skelimg=clustersedit+skelimg;
                %figure;imshow(clustersedit);
                
                endpointImage = bwmorph(clustersedit,'endpoints');
                [idx,idy]=find(clustersedit==1);
                co=[idx,idy];
                [rows, columns] = find(endpointImage);
                lr=length(rows);
                if lr>2 || lr==0 || lr==1
                    %count=count-1;
                else
                    
                    xi(1)=rows(2);yi(1)=columns(2); %endpoint as initial guess
                    xi_save=xi(1);
                    yi_save=yi(1);
                    
                    idx_find=find(co(:,1)==xi_save);
                    idy_find=find(co(:,2)==yi_save);
                    
                    [val,pos]=intersect(idx_find,idy_find);
                    
                    co_save(1,:)=[idx(val),idy(val)];
                    co(val,:)=[];
                    
                    dist_save(1)=0;
                    dist_save2(1)=0;
                    height(1)=A(idx(val),idy(val));
                    % [valy,posy]=intersect(idy_find,idx_find)  %unhash for
                    % checking purposes
                    
                    for iii = 2:length(idx)
                        [k,dist] = dsearchn(co,co_save(iii-1,:));
                        dist_save(iii)=dist+dist_save(iii-1);
                        co_save(iii,:)=co(k,:);
                        x=co(k,1);y=co(k,2);
                        height(iii)=A(x,y);
                        co(k,:)=[];

                    end
                    
                xx= smooth(co_save(:,1),prof_smoo);
                yy= smooth(co_save(:,2),prof_smoo);
                co_save_s = [xx,yy];
                 co2 = [xx,yy];         
                 co2(val,:)=[]; 
                 
                 for iii = 2:length(idx)
                        dist2 = norm(co_save_s(iii-1,:)-co_save_s(iii,:));
                        dist_save2(iii)=dist2+dist_save2(iii-1);   
                 end
                     for iii = 1:length(idx)
                        height2(iii)=A(round(co_save_s(iii,1)),round(co_save_s(iii,2)));      
                     end

                     Fco_saved{i} =co_save_s;
                     
%%width Profile analysis with filament straightening  
if wd_analysis == 1
prof=[];
                     for j = 2:numel(Fco_saved{i}(:,1))
xd = Fco_saved{i}(j,1) - Fco_saved{i}(j-1,1);
yd =Fco_saved{i}(j,2) - Fco_saved{i}(j-1,2);
    for jj=1:pleng
    perp{j}(jj,1)= Fco_saved{i}(j,1)-(jj-pleng/2)*yd;
    perp{j}(jj,2)= Fco_saved{i}(j,2)+(jj-pleng/2)*xd;
            if round(perp{j}(jj,1))>0 && round(perp{j}(jj,2))< numel(A(1,:)) && round(perp{j}(jj,2))>0 && round(perp{j}(jj,1))<numel(A(:,1)) 
                %prof(j-1,jj) = A(round(perp{j}(jj,1)),round(perp{j}(jj,2)));
                prof(j-1,jj) = F({perp{j}(jj,1),perp{j}(jj,2)});
            else
            prof(j-1,jj) = 0;
            end
            
    end
ang(j)= atan2(yd, xd);
end

 %win = 1:10:numel(Fco_saved{i}(:,1));  %width is average of 10 profiles 
%for j = 1:round(numel(Fco_saved{i}(:,1))/10)-1
%[pks(j), locs(j), wds(j)] = findpeaks(mean(prof(win(j):win(j+1),:)),'MinPeakProminence',thres,'MinPeakDistance',pleng/2);
%end 
[pks, locs, wds] = findpeaks(mean(prof(2:numel(Fco_saved{i}(:,1))-1,:)),'MinPeakProminence',thres);
%end 
                    %save prof for stright filament profile 
                    %%check that the width is measured for the correct peak
                    for w = 1:numel(wds) 
                        if locs(w) < pleng/2+pleng*loctor && locs(w) > pleng/2-pleng*loctor
                    Fwidth(i) = wds(w);  
                    Fwidth_h(i) = pks(w); 
                        else 
                    Fwidth(i) = 0; 
                    Fwidth_h(i) = 0; 
                        end
                        
                    end
else
end
                    Fpos{i}=[dist_save2]';
                    Fheight{i}=[height2]';
                    Fheightm(i) = mean(Fheight{i});
                    Flength(i)=[max(dist_save2)]';
                    dig_straight{i} = prof';
                    
      clear prof co co2 co_save co_save2 height2 height idx idy dist_save2  ang win 

                end
             end
          end
figure('Position',[10 1000 500 400])         
figure(1)

imagesc(A)
colormap(lutafm)
hold on    
for i = 1:numel(Fco_saved)
     if any(Fco_saved{i})
plot(Fco_saved{i}(:,2),Fco_saved{i}(:,1),'LineWidth',2)
str = num2str(i);
text(Fco_saved{i}(10,2),Fco_saved{i}(10,1),str,'FontSize',14,'Color','w')
     else
     end
     if exist('perp','var') == 1 && show_cross ==1
     
         for jj= 1:4:numel(perp)
         if any(perp{jj}(:))
    plot(perp{jj}(:,2),perp{jj}(:,1))
        
         end
         end
     end
end
hold off

if numel(Fco_saved)>1
    pl = 4;
else
    pl=2;
end
figure(2)
t = tiledlayout(pl,ceil(numel(Fco_saved)/2), 'Padding', 'none', 'TileSpacing', 'none'); 
for i = 1:numel(Fco_saved)
nexttile(2*i-1)
sz_dig = size(dig_straight{i});
 imagesc(dig_straight{i})
 title(i)
set(gca,'YDir','normal')
 colormap(lutafm)
 nexttile(2*i)
plot(Fheight{i})
xlim([0 numel(Fheight{i})])
end

