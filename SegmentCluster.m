classdef SegmentCluster
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        CenterPoint;
        SegmentLength;
        DirectionAngle;
        EndPoints; 
        VentLocation;
        Age;
        Error;
    end
    
    methods
        function obj = SegmentCluster(varargin)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here

            if(nargin==0)
                return;
        elseif(nargin<3)
                % input is a 2xN set of coordinates to fit a Segment
                points = varargin{1};
                
                obj.VentLocation = points;
                
                if(size(points,2)<5)
                    % not many points, treat the segment as a singleton
                    % point
                    obj.CenterPoint = mean(points,2);
                    obj.SegmentLength = 0;
                    obj.DirectionAngle = 0;
                    obj.EndPoints = nan(2,2);
                    obj.Error = obj.Distance(points);
                else
                    
                    % compute the axes
                    C = cov(points');
                    [~,s,v]=svd(C);
                    SVratio = sqrt(s(1,1)/s(2,2));
                    if(SVratio<1.5 || s(2,2)==0)
                        % points are collinear, or not very well defined
                        % as line. Treat as a point.
                        obj.CenterPoint = mean(points,2);
                        obj.SegmentLength = 0;
                        obj.DirectionAngle = 0;
                        obj.EndPoints = nan(2,2);
                    else
                        mu = mean(points,2);
                        obj.CenterPoint = mu;
                        obj.DirectionAngle = atan2(v(2,1),v(1,1));
                        
                        DistFromCenter = sqrt(sum(bsxfun(@minus,points,mu).^2));
                        obj.SegmentLength = 3*mean(DistFromCenter);
                        
                        % enforce maximum length
                        if(obj.SegmentLength>10000)
                            obj.SegmentLength = 10000;
                        end
                        
                        obj.EndPoints = [mu mu]+obj.SegmentLength/2*...
                            [-v(:,1), v(:,1)];
                    end
                    
                    obj.Error = obj.Distance(points);
                end
                
                if(nargin==2)
                    obj.Age = mean(varargin{2});
                end
            else
                % three inputs that define the properties
                obj.CenterPoint = varargin{1};
                obj.SegmentLength = varargin{2};
                obj.DirectionAngle = varargin{3};
                obj.VentLocation = varargin{4};
                obj.Age = varargin{5};
                v = [cos(obj.DirectionAngle); sin(obj.DirectionAngle)];
                obj.EndPoints = [obj.CenterPoint, obj.CenterPoint]+obj.SegmentLength/2*...
                    [-v(:,1) v(:,1)];
                obj.Error = obj.Distance(obj.VentLocation);
            end
        end
        
        function plotSegment(obj,varargin)
            if(obj.SegmentLength==0)
                if(nargin==1)
                    plot(obj.CenterPoint(1),obj.CenterPoint(2),'x');
                else
                    plot(obj.CenterPoint(1),obj.CenterPoint(2),varargin{1});
                end
            else
                
                dir = [cos(obj.DirectionAngle); sin(obj.DirectionAngle)];
                X = obj.CenterPoint+obj.SegmentLength/2*[-dir, dir];
                if(nargin==1)
                    plot(X(1,:),X(2,:));
                else
                    plot(X(1,:),X(2,:),varargin{1});
                end
            end
        end
        
        function D = Distance(obj,Points)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if(obj.SegmentLength==0)
                % Segment Cluster defined as a point
                D = sqrt(sum(bsxfun(@minus,Points,obj.CenterPoint).^2));
            else
                R = [cos(obj.DirectionAngle), -sin(obj.DirectionAngle); ...
                    sin(obj.DirectionAngle) cos(obj.DirectionAngle)];
                X = R'*bsxfun(@minus,Points,obj.CenterPoint);
                
                D = zeros(1,size(Points,2));
                left = find(X(1,:)<=-obj.SegmentLength/2);
                if(~isempty(left))
                    D(left) = sqrt(X(2,left).^2+...
                        (X(1,left)+obj.SegmentLength/2).^2);
                end
                
                right = find(X(1,:)>=obj.SegmentLength/2);
                if(~isempty(right))
                    D(right) = sqrt(X(2,right).^2+...
                        (X(1,right)-obj.SegmentLength/2).^2);
                end
                
                center = find(X(1,:)>-obj.SegmentLength/2 & ...
                    X(1,:)<obj.SegmentLength/2);
                if(~isempty(center))
                    D(center) = abs(X(2,center));
                end
            end
        end
    end
end