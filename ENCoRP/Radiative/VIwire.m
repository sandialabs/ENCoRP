function Vf = VIwire(xf,yf,zf,wpwf,extraf,ff,VorCf)

    % x = x (m)
    % y = y (m)
    % z = z (m)
    % wwv = whichphasewire
    % extra = 1, 0, -1 = forward derivative, no derivative, backward derivative
    % f = frequency
    % VorC = voltage 1 or current 2
    
    h = 0.001; % For derivative

    deriv = false;
    if extraf(1) ~= 0
        deriv = true;
    end
    
    %
    % Load wires and nodal locations
    %
    load('GUIData.mat','GUIData');
    Freq = GUIData.Freq;
    VI = GUIData.VI;
    ncond = size(VI{1,1},2)/2;
   
    %
    % Pull useful variables
    %
    xv = GUIData.Nodes(:,1)';
    yv = GUIData.Nodes(:,2)';
    zv = GUIData.Nodes(:,3)';
    W = GUIData.WireID;

    Vf = zeros(length(xf),1);
    for i2 = 1:length(xf)
        x = xf(i2);
        y = yf(i2);
        z = zf(i2);
        f = ff(i2);
        wpw = wpwf(i2);
        extra = extraf(i2);
        VorC = VorCf(i2);

        % Acquire frequency ind   
        [~,nind] = min(abs(f-Freq));

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
        % Now pull the voltages/currents corresponding to the wire (modal/mean)
        %
        freqv = nind; 
        if VorC == 1
            Voltr = real(squeeze(VI{freqv,wind}(:,wpw)));
            Volti = imag(squeeze(VI{freqv,wind}(:,wpw)));
        elseif VorC == 2
            Voltr = real(squeeze(VI{freqv,wind}(:,wpw + ncond)));
            Volti = imag(squeeze(VI{freqv,wind}(:,wpw + ncond)));
        else
            error('Nice coding... VorC.')
        end
    
        %
        % Translate to one dimension and interpolate
        %
        xvi = linspace(0,1,length(Voltr))';
        Pvi = norm(P' - [xv(W(wind,1)),yv(W(wind,1)),zv(W(wind,1))]) / norm([xv(W(wind,2)),yv(W(wind,2)),zv(W(wind,2))] - [xv(W(wind,1)),yv(W(wind,1)),zv(W(wind,1))]);
        Vr = interp1(xvi,Voltr,Pvi);
        Vi = interp1(xvi,Volti,Pvi);
        V = Vr+1j*Vi;
        
        %
        % Compute derivative if requested
        %
        if deriv
            % Normalize h by wire length
            h = h/norm([xv(W(wind,2)),yv(W(wind,2)),zv(W(wind,2))] - [xv(W(wind,1)),yv(W(wind,1)),zv(W(wind,1))]);
            if extra == 1
                Vr = interp1(xvi,Voltr,min(1,Pvi+h));
                Vi = interp1(xvi,Volti,min(1,Pvi+h));
            elseif extra == -1
                Vr = interp1(xvi,Voltr,max(0,Pvi-h));
                Vi = interp1(xvi,Volti,max(0,Pvi-h));
            end
            if extra == 1 || extra == -1 % last chance to avoid
                V = (Vr+1j*Vi).*(([xv(W(wind,2)),yv(W(wind,2)),zv(W(wind,2))] - [xv(W(wind,1)),yv(W(wind,1)),zv(W(wind,1))])/norm([xv(W(wind,2)),yv(W(wind,2)),zv(W(wind,2))] - [xv(W(wind,1)),yv(W(wind,1)),zv(W(wind,1))]));
            end
        end
        if extra ~= 0
            Vf = V;
            return
        else
            Vf(i2) = V;
        end
    end
    Vf(isnan(Vf)) = 0;
    Vf(isinf(Vf)) = 0;
end
