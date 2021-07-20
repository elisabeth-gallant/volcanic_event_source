function Z = linkage_Segment(d)
    M = size(d,2);
    Z = inf(M-1,3);
    
    % initial line segment objects
    for m=1:M
        ActiveClusters(m) = SegmentCluster(d(:,m));
    end
    NActive = M;
    CurrentIndex = 1:M;
    for m=1:M-1
        Dmat = inf(NActive);
        
        % construct distance matrix
        for k=1:NActive
            for j=k+1:NActive
                Dmat(k,j) = ...
                    Segment_Cluster_Distance(ActiveClusters(k),ActiveClusters(j));
            end
        end
        
        MergeDistance = min(Dmat(:));
        [IX, IY] = find(Dmat==MergeDistance,1,'first');
        if(CurrentIndex(IX)<CurrentIndex(IY))
            Z(m,:) = [CurrentIndex(IX), CurrentIndex(IY), MergeDistance];
        else
            Z(m,:) = [CurrentIndex(IY), CurrentIndex(IX), MergeDistance];
        end
        ActiveClusters(end+1) = SegmentCluster([ActiveClusters(IX).VentLocation, ...
            ActiveClusters(IX).VentLocation]);
        CurrentIndex(end+1) = m+M;
        CurrentIndex([IX, IY]) = [];
        ActiveClusters([IX,IY]) = [];
        
        NActive = NActive-1;
    end
end