%This code was last updated on 19 July, 2021 by Lis Gallant

% you will need to change the following lines to use your own data:
% 11 = your background DEM file for plotting purposes)
% 24 = your excel file with your vent locations and ages. If you have a text
%  file instead of an excel file, comment out line 21 and uncomment 24, 25
% 64 = change those values to suit your specific volcanic field
% 115 = change output file name for event coordinates to be saved in .csv

%This line provides the background map image for the plots
DMFile = '../DEMs/idaho_map_utm_hillshade.tif';

DEM = readgeoraster(DMFile);
DEMinfo = geotiffinfo(DMFile);

IMG = flipud(DEM);
IMG(IMG(:)==-9999) = min(IMG(IMG(:)~=-9999));

X = linspace(DEMinfo.BoundingBox(1,1),DEMinfo.BoundingBox(2,1),size(DEM,1));
Y = linspace(DEMinfo.BoundingBox(1,2),DEMinfo.BoundingBox(2,2),size(DEM,2));

%vents input: this reads an excel file and selects data based on tab name
vents = xlsread('../vent_data_all.xls','Craters of the Moon');

%vents = load('../Pali-Aike_All.txt');
%vents(:,3) = 1;

Nv = size(vents,1);
[AllVentGMM, BaselineBW] = GaussKDE(vents(:,1:2)');

AX = [min(vents(:,1:2)); max(vents(:,1:2))];
AX(1,:) = AX(1,:)-20000;
AX(2,:) = AX(2,:)+20000;

% no clustering
Ns = 1000;
Xs = linspace(AX(1,1),AX(2,1),Ns)';
Ys = linspace(AX(1,2),AX(2,2),Ns)';
PDFALL = zeros(Ns,Ns);
for k=1:Ns
    PDFALL(k,:) = AllVentGMM.pdf([Xs,Ys(k)*ones(Ns,1)]);
end
image(X,Y,repmat(double(IMG)/double(max(IMG(:))),[1,1,3]));
set(gca,'ydir','normal'); hold on;
contour(Xs,Ys,PDFALL);

%% segment based clustering 
% This is the line that you want to change to get values specific to your
% volcanic field: [Events, ID] = LSBE(vents, #1, #2, min(#3, size(vents,1)
% #1 corresponds to the temporal cut-off value - I choose a value closer to
% zero when I have higher confidence in the data, such as with Yucca
% Mountain
% #2 corresponds to the maximum distance between your vents within the
% clusters - this is specific to each volcanic field based on the tectonic
% setting, but somewhere between 5000 and 10000 meters is a good start
% #3 corresponds to the maximum number of clusters you're willing to
% consider - for a volcanic field like Pali-Aike, with a high number of
% vents, this number is important because it gives you a bit more control
% over how much lumping or splitting the code does. If you'd like to
% consider all available outputs, set this value to the # of vents within your
% field


[Events, ID] = Line_Segment_Based_Events(vents,0.7,2500,min(72,size(vents,1)));

% Yucca Mountain inputs: [Events, ID] = Line_Segment_Based_Events(vents,0.5,4000,min(42,size(vents,1)));
% Craters of the Moon inputs: [Events, ID] = Line_Segment_Based_Events(vents,0.7,2500,min(72,size(vents,1)));
% Pali Aike inputs: [Events, ID] = Line_Segment_Based_Events(vents,0.5,5000,min(450,size(vents,1)));
% The higher the number in for the maximum number of clusters (#3), the
% longer the run time.

MinVents = 3;
Mu = zeros(Nv,2);
Sig = zeros(2,2,Nv);
Ne = length(Events);
Mix = ones(Nv,1)/Ne;
n=1;
PDFEvent = zeros(Ns,Ns);
for k=1:Ne
    CurNumVents = size(Events(k).VentLocation,2);
    ix = n:(n+CurNumVents-1);
    if(CurNumVents<MinVents)
        % not enough vents in cluster to want to compute cluster pdf, use
        % the baseline
        Mu = Events(k).VentLocation';
        Sig = repmat(BaselineBW,[1,1,CurNumVents]);
        EventGMM = gmdistribution(Mu,Sig);
    else
        [EventGMM, BWEvent] = GaussKDE(Events(k).VentLocation);
    end
   
    for m=1:Ns
        PDFEvent(m,:) = PDFEvent(m,:)+EventGMM.pdf([Xs,Ys(m)*ones(Ns,1)])';
    end
end
PDFEvent = PDFEvent/Ne;

figure;
image(X,Y,repmat(double(IMG)/double(max(IMG(:))),[1,1,3]));
set(gca,'ydir','normal'); hold on;
plot(vents(:,1),vents(:,2),'.r')
% axis([287000 317000 4760000 4840000])
for k=1:Ne
    if(~isempty(Events(k).EndPoints))
    plot(Events(k).EndPoints(1,:),Events(k).EndPoints(2,:),'b')
    plot(Events(k).CenterPoint(1),Events(k).CenterPoint(2),'ob')
    end
end

%This section pumps out the X,Y coordinates of each events
EventCenter = zeros(length(Events),2);
for k=1:length(Events)
EventCenter(k,:) = Events(k).CenterPoint;
end
 writematrix(EventCenter, 'com.csv');
 function EventTable = write_event_csv(Events,filename)

N = length(Events);

Array = [];
for n=1:N
    Nv = size(Events(n).VentLocation,2);
    for k=1:Nv
        Array = [Array; Events(n).VentLocation(:,k)', n, Events(n).CenterPoint', Events(n).Age];
    end
end

varnames = {'Vent_N', 'Vent_E', 'Event_ID', 'Event_N', 'Event_E', 'Event_Age'};

EventTable = array2table(Array,'VariableNames',varnames);

writetable(EventTable,filename);
 end