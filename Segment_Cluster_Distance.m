function distance = Segment_Cluster_Distance(C1,C2)
    
    % flags indicating if Cluster is defined by a point
    Point1 = any(isnan(C1.EndPoints(:)));
    Point2 = any(isnan(C2.EndPoints(:)));
    
    if(Point1 && Point2)
        % both are points
        distance = norm(C1.CenterPoint-C2.CenterPoint);
    elseif((~Point1) && (~Point2))
        % both are line segments
        t = linspace(0,1,100);
        NT = length(t);
        F = zeros(1,NT);
        dL = C1.EndPoints(:,2)-C1.EndPoints(:,1);
        Line = C2.EndPoints;
        for m=1:NT
            pt = C1.EndPoints(:,1)+t(m)*dL;
            for n=1:NT
                f = sqrt((pt(1)-((Line(1,2)-Line(1,1))*t+Line(1,1))).^2+...
                    (pt(2)-((Line(2,2)-Line(2,1))*t+Line(2,1))).^2);
            end
            F(m) = trapz(t,f);
        end
        distance = trapz(t,F);
    else
        % one is a point the other a line segment
        if(Point1)
            pt = C1.CenterPoint;
            Line = C2.EndPoints;
        else
            pt = C2.CenterPoint;
            Line = C1.EndPoints;
        end
        % numerically integrate to find the distance
        t = linspace(0,1,100);
        F = sqrt((pt(1)-((Line(1,2)-Line(1,1))*t+Line(1,1))).^2+...
            (pt(2)-((Line(2,2)-Line(2,1))*t+Line(2,1))).^2);
        distance = trapz(t,F);
    end
    if(distance==0)
        distance = distance;
    end
end