% vents(N,3) - East, North, Time (NaN indicates no time information)
% template(3,1) - major axis half length, minor axis half length,
%   orentation angle (deg)
% TimeCutoff - clustering threshold for time clustering
function [Events, ID] = Line_Segment_Based_Events(vents,TimeCutoff,minDist,maxClust)
    Events = SegmentCluster.empty;
    
    ID = zeros(size(vents,1),1);
    
    % extract vents with ages
    ixAgeVent = ~isnan(vents(:,3));
    ENTAgeAll = vents(ixAgeVent,:);
    ENTNoAgeAll = vents(~ixAgeVent,:);
    
    % perform the clustering of the Ages
    [ENTime, ~, Tc] = Time_Clustering(ENTAgeAll,TimeCutoff);
    
    % perform spatial clustering via the template for each time cluster
    total_class = 0;
    AgeAssignment = zeros(sum(ixAgeVent),1);
    
    for k=1:max(Tc)
        % indices to vents in the time cluster
        inTimeCluster = find(Tc==k);
        
        % vents in the time cluster
        VentsTimeCluster = ENTime{k}.vents(:,1:2);
        
        % check if the time clustering has more than 1 data point
        if(size(VentsTimeCluster,1)>1)
            mu = mean(VentsTimeCluster);
        else
            mu = VentsTimeCluster;
        end
        Dist = sqrt(sum((VentsTimeCluster-repmat(mu,size(VentsTimeCluster,1),1)).^2,2));
        if(all(Dist<=minDist))
            seg = SegmentCluster();
            seg.Error = Dist;
            seg.VentLocation = VentsTimeCluster';
            seg.DirectionAngle = 0;
            seg.SegmentLength = 0;
            seg.EndPoints = nan(2);
            seg.CenterPoint = mu';
            seg.Age = mean(ENTime{k}.vents(:,3));
            
            % points are packed tight, only 1 cluster
            Ts = ones(length(inTimeCluster),1);
            
        else
            M = min(maxClust,ceil(length(inTimeCluster)/2)-1);
            if(M<2)
                M=2;
            end
            % points have enough spread to consider multiple clusters
            [seg, Ts, ~, DBI] = Hierarchical_Segment(VentsTimeCluster',ENTime{k}.vents(:,3)',M);
        end
        
        % Add results to output (nspace = # of clusters)
        AgeAssignment(inTimeCluster) = Ts+total_class;
        total_class = total_class+max(Ts);
        Events = [Events, seg];
    end
    
    ID(ixAgeVent) = AgeAssignment;
    
    % cluster the vents without ages
    N = sum(~ixAgeVent);
    
    
    if(N>1)
        mu = mean(ENTNoAgeAll(:,1:2));
    else
        mu = ENTNoAgeAll(:,1:2);
    end
    Dist = sqrt(sum((ENTNoAgeAll(:,1:2)-repmat(mu,N,1)).^2,2));
    
    if(N==1 || all(Dist<=minDist))
        % only 1 vent without an age
        total_class = total_class+1;
        NoAgeAssignment = total_class*ones(N,1);
        % Add results to output
        seg = SegmentCluster();
        seg.Error = Dist;
        seg.VentLocation = ENTNoAgeAll(:,1:2)';
        seg.DirectionAngle = 0;
        seg.SegmentLength = 0;
        seg.EndPoints = nan(2);
        seg.CenterPoint = mu';
        seg.Age = mean(ENTNoAgeAll(:,3));
        
        Events = [Events, seg];
    elseif(N~=0)
        M = min(maxClust,N-1);
        % points have enough spread to consider multiple clusters
        [seg, Ts] = Hierarchical_Segment(ENTNoAgeAll(:,1:2)',ENTNoAgeAll(:,3)',M);
        
        NoAgeAssignment = Ts+total_class;
        total_class = total_class+max(Ts);
        
        Events = [Events, seg];
    else
        % no ageless vents
        NoAgeAssignment = [];
    end
    
    ID(~ixAgeVent) = NoAgeAssignment;
    

ix = false(1,length(Events));
for n=1:length(Events)
    if(isempty(Events(n).CenterPoint))
        ix(n) = true;
    end
end
Events(ix) = [];

end