function [ENTime, Z, Tc] = Time_Clustering(ENAge,cutoff)
    
    T = ENAge(:,3);
    
    % construct a heirarchical tree describing the time data
    Z = linkage(T);
    
    % perform temporal clustering one vents with ages
    Tc = cluster(Z,'cutoff',cutoff);
    
    NC = max(Tc);
    ENTime = cell(NC,1);
    % calculate the mean of each cluster and store the members of each cluster
    for k=1:NC
        curr = Tc==k;
        ENTime{k}.time = mean(T(curr));
        ENTime{k}.vents = ENAge(curr,:);
        ENTime{k}.NData = sum(curr);
    end
end