function [seg, selClust, Error, DBIset, LINK] = Hierarchical_Segment(z,Time,maxClust)
    
    DBIset = inf(1,maxClust);
    LINK = linkage_Segment(z);
    
    for n=2:maxClust
        % cluster assuming n clusters
        Tc = cluster(LINK,'maxclust',n);
        
        % compute DBI
        for k=1:n
            ActiveClusters(k) = SegmentCluster(z(:,Tc==k));
        end
        
        DBI = zeros(n,n);
        for k=1:n
            for j=k+1:n
                intraCluster = Segment_Cluster_Distance(ActiveClusters(k),ActiveClusters(j));
                DBI(k,j) = (mean(ActiveClusters(k).Error)+mean(ActiveClusters(j).Error))/intraCluster;
            end
        end
        DBI = DBI+DBI';
        DBIset(n) = mean(max(DBI));
    end
    
    [~, n] = min(DBIset);
    selClust = cluster(LINK,'maxclust',n);
    Error = 0;
    for k=1:n
        seg(k) = SegmentCluster(z(:,selClust==k),Time(selClust==k));
        Error = Error+sum(seg(k).Error);
    end
end