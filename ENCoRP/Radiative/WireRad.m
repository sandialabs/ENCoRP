function rw = WireRad(xf,yf,zf,wpwf)

    % x = x (m)
    % y = y (m)
    % z = z (m)
    % wwv = whichphasewire

    %
    % Load wires and nodal locations
    %
    load('GUIData.mat','GUIData');
   
    %
    % Pull useful variables
    %
    xv = GUIData.Nodes(:,1)';
    yv = GUIData.Nodes(:,2)';
    zv = GUIData.Nodes(:,3)';
    W = GUIData.WireID;
    OldW = GUIData.WiresOld;
    rad = GUIData.Radii;

    rw = zeros(length(xf),1);
    for i2 = 1:length(xf)
        x = xf(i2);
        y = yf(i2);
        z = zf(i2);
        wpw = wpwf(i2);

        %
        % Find nearest wire
        %
        P = [x,y,z]';
        ww = -ones(size(W,1),1);
        for i = 1:size(W,1)
            V1 = [xv(W(i,1)),yv(W(i,1)),zv(W(i,1))]';
            V2 = [xv(W(i,2)),yv(W(i,2)),zv(W(i,2))]';
            d1 = norm(P-V1);
            d2 = norm(P-V2);
            %
            % Compute shortest distance to line
            %
            d3 = V1-V2;
            d4 = P-V2;
            if any(size(d3)~=size(d4))
                sz1 = size(d3);
                sz2 = size(d4);
                str = sprintf('%d,%d,%d,%d',sz1(1),sz1(2),sz2(1),sz2(2));
                error(str); %#ok<*SPERR> 
            end
            dl = norm(cross(d3,d4))/norm(d3);
            %
            % Compare to line segment not infinite line
            %
            if max(d1,d2) <= norm(V1-V2)+1e-1
                ww(i) = dl;
            else
                ww(i) = Inf;
            end
        end
    
        %
        % Now choose nearest wire
        %
        [checkval,wind] = min(ww);
        tol = 1e-1;
        if checkval > tol
            error(sprintf('%f,%f,%f,%f',checkval,P(1),P(2),P(3)))
        end
        
        %
        % Match nearest new wire to old wire
        %
        for i3 = 1:size(OldW,1)
            if OldW(i3,1) == W(wind,1) && OldW(i3,2) == W(wind,2)
                rw(i2) = rad(i3,wpw+1);
                break
            end
            if OldW(i3,2) == W(wind,1) && OldW(i3,1) == W(wind,2)
                rw(i2) = rad(i3,wpw+1);
                break
            end
        end
    end
    rw(isnan(rw)) = 0.001;
    rw(isinf(rw)) = 0.001;
    if any(rw==0)
        error('Zero rw.')
    end
end
