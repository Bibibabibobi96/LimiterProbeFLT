function [in_or_out]=in_or_out_area(divertor_r,divertor_z,r,z)
% clear
% clc
% the sequence points of an area. the last point is the same as the first one
if isnan(r) | isnan(z)
    in_or_out = 0;
else
    areapoints(:,1) = divertor_r;
    areapoints(:,2) = divertor_z;
    
    % plot(areapoints(1:end,1),areapoints(1:end,2));
    n_area=length(areapoints)-1;
    lines=[areapoints(1:end-1,:) areapoints(2:end,:)]; % all the border lines. save the two points of a line
    
    % lines
    % for: the longitude of the first point must be less or equal to the second.
    for n_point=1:n_area
        if lines(n_point,1)>lines(n_point,3)
            lines(n_point,:)=[lines(n_point,3:4),lines(n_point,1:2)]; % change the location of the two points of a line.
        end
    end
    
    % lines
    % hold on;
    ma=max(areapoints); % max value of longitude and latitude of the area.
    mi=min(areapoints); % min value of longitude andlatitude of the area.
    % testp=mi+(ma-mi).*rand(1,2); % test point.
    
    testp=[r,z];
    %testp=[116.293859, 40.983803];
    % plot(testp(1),testp(2),'*r'); % plot the testpoint
    count=0; % the number of intersective points.
    % if: if the point is out of the min square of the area, the test point must be
    % out of the area and we can return directly.
    if testp(1)>ma(1) || testp(1)<mi(1) || testp(2)>ma(2) || testp(2)<mi(2)
        %     disp('out_0');
        in_or_out = 0;
        return;
    end
    
    % for: compute the count of every line.
    for i=1:n_area
        if testp(1)>lines(i,1) && testp(1)<lines(i,3) % the longitude of the point is between the two points of the line.
            y=lines(i,2)+(lines(i,4)-lines(i,2))/(lines(i,3)-lines(i,1))*(testp(1)-lines(i,1));
            if y<testp(2) % have one intersective point with the line.
                count=count+1;
            end
            if y==testp(2)  % on the border of the area. and return.
                %         disp('in');
                in_or_out = 1;
                return;
            end
        end
    end
    
    % for: compute the count of every vetex. only the left side (the right side is also ok).
    for i=1:n_area
        if lines(i,1)==testp(1) && lines(i,2)<testp(2) && lines(i,3)>test(1) % have one intersective point with the vetex. and do not consider the vetrical line of the area.
            count=count+1;
        elseif lines(i,1)==testp(1) && lines(i,2)==testp(2) % on the vetex.
            %         disp('in');
            in_or_out = 1;
        end
    end
    
    if mod(count,2)==1
        %     disp('in');
        in_or_out = 1;
    else
        %     disp('out');
        in_or_out = 0;
    end
end
end