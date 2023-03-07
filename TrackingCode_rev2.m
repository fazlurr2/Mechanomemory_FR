
%%% Fazlur Rashid, August 2021

clc;
clear;
close all;


global Counter
pathname='C:\';
pixelsize=2;
threshold=2; 

mag = 0.103*1000; %%in nm   %40X objective: 0.163 um/pixel; 63X objective: 0.103 um/pixel


[filename1,pathname]=uigetfile([pathname,'*.*'],'select first image in sequence');
[filename2,pathname]=uigetfile([pathname,'*.*'],'select last image in sequence');


cd(pathname);
mkdir(pathname,'Analysis');

I=imread([pathname,filename1]);

[I2, rect]=imcrop(I);
imshow(I2);

k1 = strfind(filename1, '_');
k2 = strfind(filename1, '+');
k3= strfind(filename1, '-');
k=[k1 k2 k3];
k=sort(k);
k1=k(1);
k2=k(2);


CounterStart=str2num(filename1([k1+1:k2-1]));
CounterEnd=str2num(filename2([k1+1:k2-1]));
Counter=0;


PositionList=[];
mkdir([pathname,'Analysis'],'Crop');
mkdir([pathname,'Analysis'],'Centroid');

cd([pathname,'Analysis']);
writerObj = VideoWriter('Video.avi');
open(writerObj);


% % % compensate for drift between first and last image
DriftComp=1;   %Turn Drift Comenseation ON (==1) or OFF(==0)
dxy=[0 0];
if DriftComp==1

frames=0;
for i=[CounterStart:CounterEnd]
    filename=['*',num2str(i),'*.*'];
    cd(pathname);
    filenameX=dir([filename]);          %% Find images in range
    
    if isempty(filenameX)
        continue
    end
    
    frames=frames+1;
end

clf;

% First Frame
I=imread([pathname,filename1]);
I2=imcrop(I,rect);
I2b=bpass(I2,1,pixelsize+1);


pk = pkfnd(I2b,threshold,pixelsize); 
cnt1 = cntrd(I2b,pk,pixelsize+2);
figure(1);subplot(1,2,1)
colormap('gray'); imagesc(I2b);axis square;

for i=1:size(cnt1,1)
    hold on; text(cnt1(i,1),cnt1(i,2)-3,num2str(i),'Fontsize', 16,'Color', 'red')
end



I=imread([pathname,filename2]);
I2=imcrop(I,rect);
I2b=bpass(I2,1,pixelsize+1);
pk = pkfnd(I2b,threshold,pixelsize); 
cnt2 = cntrd(I2b,pk,pixelsize+2);
subplot(1,2,2); colormap('gray'); imagesc(I2b);axis square;



for i=1:size(cnt2,1)
    hold on; text(cnt2(i,1),cnt2(i,2)-3,num2str(i),'Fontsize', 16,'Color', 'blue')
end
figure(1); 
im1=input('enter the particle number for drift from image 1 (in array format)');

im2=input('enter the particle number for drift from image 2 (in array format)');

if (length(im1)==1) || (length(im2)==1)
    dxy=((cnt2(im2,[1,2])-cnt1(im1,[1,2]))./frames);
else
    dxy=mean((cnt2(im2,[1,2])-cnt1(im1,[1,2]))./frames);   %Average drift per frame
end

end
%%%%---------------end drift compensation


time=[];

for i=[CounterStart:CounterEnd]
    filename=['*',num2str(i),'*.*'];
    cd(pathname);
    filenameX=dir([filename]);          %% Find images in range
    
    if isempty(filenameX)
        continue
    end
    
    time=[time i];
    Counter=Counter+1;
    filename=filenameX.name;
    I=imread([pathname,filename]);
    I2=imcrop(I,rect);
    
    cd([pathname,'Analysis','/Crop']);
    imwrite(I2,[filename1(1:k1),num2str(i),'_crop.tiff'],'tiff')
    
%%%%Filter and find centroids of particles
    I2b=bpass(I2,1,pixelsize+1);
    figure(i); colormap('gray'); imagesc(I2b)
    
    pk = pkfnd(I2b,threshold,pixelsize); 

    cnt = cntrd(I2b,pk,pixelsize+2);
% % % %     cnt=pk;
    
    if isempty(cnt)
        continue
    end

    drift=(ones(size(cnt,1),1)*dxy*Counter);
    cnt(:,[1,2])=cnt(:,[1,2])-drift;   %%%compensated for drift


    figure(i);hold on;plot(cnt(:,1),cnt(:,2),'ro');
    cd([pathname,'Analysis','/Centroid']);
    saveas(gcf,[filename1(1:k1),num2str(i),'_Centroid.fig']);
    frame = getframe;
    writeVideo(writerObj,frame);
    close(gcf)

    
    s=size(cnt);
    PositionList=[PositionList;[cnt(:,[1,2]) Counter*ones(s(1),2)]];
end
close(writerObj);

%%% Find the spots with continuous trajectories
  
tr = track(PositionList,3); %set paprameters in the code track.m 
%%param.good set such that only tracks spots that exist in the whole
%%imaging sequence

tr(:,[1,2])=tr(:,[1,2])*mag;

cd([pathname,'Analysis\Centroid']);
open([filename1(1:k1),num2str(CounterStart),'_Centroid.fig'])
hold on;

TrackNum=tr(end,end);
TrackNum1=1;
% TrackNum=1;
x0=0;
y0=0;

disp=[];
for i=TrackNum1:TrackNum
    row=find(tr(:,4)==i);
    x=tr(row,1)-x0;
    y=tr(row,2)-y0;
    plot(x/mag,y/mag,'.-')
    text(x(1)/mag,y(1)/mag-2,num2str(i),'Fontsize', 16,'Color', 'red')
    hold on;
    ddisp(:,i)=sqrt((x(2:end)-x(1:end-1)).^2+(y(2:end)-y(1:end-1)).^2);
    disp=[disp;[ddisp(:,i)]];
end
cd([pathname,'Analysis']);
saveas(gcf,['Trajectory.fig']);
saveas(gcf,['Trajectory.tiff']);




%%%%-------------
%%Calculate the histogram

[NN,XX]=hist(disp(:),[0:10:200]);
NN=NN/sum(NN)*100;    % percentage of displacements for given displacement
figure; bar(XX,NN,'histc')
xlabel 'displacement (nm)'
ylabel '%'
ylim([0 60])
xlim([0 200])
cd([pathname,'Analysis']);
saveas(gcf,['DisplacementHistogram.fig']);
saveas(gcf,['DisplacementHistogram.tiff']);

%%%%-------------
%%Calculate MSD
for i=TrackNum1:TrackNum
    row=find(tr(:,4)==i);
    x=tr(row,1)-x0;
    y=tr(row,2)-y0;
    
    for tau=[1:length(x)-1]
        msdd=[];
        for t=1:length(x)-tau
            msdd=[msdd (x(t+tau)-x(t)).^2+(y(t+tau)-y(t)).^2];
        end
        MSD(tau,i)=mean(msdd);
    end   
      
end

time=(time-time(1))/10;
time=time(:, [1:size(MSD,1)+1]);

figure;plot(time(1:end-1),MSD)
xlabel 'time (s)'
ylabel 'MSD (nm^2)'
cd([pathname,'Analysis']);
saveas(gcf,['MSD.fig']);
saveas(gcf,['MSD.tiff']);


figure;plot(time(1:end-1),ddisp)
xlabel 'time (s)'
ylabel 'displacement (nm)'
cd([pathname,'Analysis']);
saveas(gcf,['Displacement.fig']);
saveas(gcf,['Displacement.tiff']);


pathname

cd([pathname,'Analysis']);
save AllVariable.mat

