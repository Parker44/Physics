function [t, r] = mutual_orbit(tmax, level, m1, m2, d)

    % Discretize continuum domain
    nt = 2^level + 1;
    t = linspace(0.0, tmax, nt);
    deltat = t(2) - t(1);
    
    m = m1 + m2;
    r = zeros(2, 3, nt); 
    v = zeros(2, 3, nt);
    
    r1_0 = (m2/m)*d;
    r2_0 = (m1/m)*d;
    
    r(1, :, 1) = [r1_0 0 0];
    r(2, :, 1) = [-r2_0 0 0];
    
    v(1, :, 1) = [0 sqrt(m2*r1_0)/d 0];
    v(2, :, 1) = [0 -sqrt(m1*r2_0)/d 0];
    
    a(1, :, 1) = [-m2/(d^2) 0 0];
    a(2, :, 1) = [m1/(d^2) 0 0];
    
    r(1, :, 2) = r(1, :, 1) + v(1, :, 1)*deltat + 0.5 * deltat^2 * a(1, :, 1);
    r(2, :, 2) = r(2, :, 1) + v(2, :, 1)*deltat + 0.5 * deltat^2 * a(2, :, 1);
    
    %-----------------------------------------------------------
    % Set plotenable to non-zero/zero to enable/disable plotting.
    %-----------------------------------------------------------
    plotenable = 1;
    plt = gobjects(1, 2);
    %-----------------------------------------------------------
    % Parameter to control speed of animation.  Script execution
    % will pause for pausesecs each time a new frame is drawn.
    % 
    % Setting this parameter to a "largish" value, say 0.1
    % (seconds), will produce a slow-motion effect.
    % 
    % Set it to 0 for maximum animation speed.
    %-----------------------------------------------------------
    pausesecs = 0.01;

    %-----------------------------------------------------------
    % Plot attributes defining the appearance  of the ball.
    %-----------------------------------------------------------

    % Particles have a (marker) size of 5
    bodysize = 2;
    % they are yellow
    bodycolor = 'y';
    % ... and it's plotted as a circle.
    bodymarker = 'o';

    % Plot background is black
    bgcolor = 'k';

    %-----------------------------------------------------------
    % Set avienable to a non-zero value to make an AVI movie.
    %-----------------------------------------------------------
    avienable = 1;

    % If plotting is disabled, ensure that AVI generation
    % is as well
    if ~plotenable
       avienable = 0;
    end

    % Name of avi file.
    avifilename = 'mutual_orbit.avi';

    % Presumed AVI playback rate in frames per second.
    aviframerate = 25;

    %-----------------------------------------------------------
    
    %-----------------------------------------------------------
    % If AVI creation is enabled, then initialize an avi object.
    %-----------------------------------------------------------
    if avienable
       aviobj = VideoWriter(avifilename);
       open(aviobj);
    end
    
    
    if plotenable
        % Clear figure
        clf;

        % Don't erase figure after each plot command.
        hold on;

        % Define plotting area, using a 1:1 aspect ratio for the 
        % plotted region, boxed axes and a 15%-width "border" around 
        % the unit square.
        axis square;
        xlim([-1, 1]);
        ylim([-1, 1]);
        set(gca, 'Color', bgcolor);
        set(gcf, 'Visible', 'off');
        set(gcf,'Position',[100 100 250 250])
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        
        % Draw the particles. 
        plt(1) = plot(r(1, 1, 1), r(1, 2, 1), 'Marker', bodymarker, 'MarkerSize', bodysize, ...
         'MarkerEdgeColor', bodycolor, 'MarkerFaceColor', bodycolor);

        plt(2) = plot(r(2, 1, 1), r(2, 2, 1), 'Marker', bodymarker, 'MarkerSize', bodysize, ...
         'MarkerEdgeColor', bodycolor, 'MarkerFaceColor', bodycolor);
    end
    
    for n=2:nt-1
        
        if plotenable

            % Update the particle positions. 
            plt(1).XData = r(1, 1, n);
            plt(1).YData = r(1, 2, n);
            
            plt(2).XData = r(2, 1, n);
            plt(2).YData = r(2, 2, n);

            % Record video frame if AVI recording is enabled and record 
            % multiple copies of the first frame.  Here we record 5 seconds
            % worth which will allow the viewer a bit of time to process 
            % the initial scene before the animation starts.
            if avienable
                if t == 0
                    framecount = 5 * aviframerate ;
                else
                    framecount = 1;
                end
                for iframe = 1 : framecount
                    writeVideo(aviobj, getframe(gcf));
                end
            end

            % Pause execution to control interactive visualization speed.
            pause(pausesecs);
        end

        a = twobodyaccn(m1, m2, r(:, :, n));

        r(1, :, n+1) = 2*r(1, :, n) - r(1, :, n-1) + deltat^2 * a(1, :);
        r(2, :, n+1) = 2*r(2, :, n) - r(2, :, n-1) + deltat^2 * a(2, :);

    end
    
    if plotenable
        plt(1).XData = r(1, 1, nt);
        plt(1).YData = r(1, 2, nt);

        plt(2).XData = r(2, 1, nt);
        plt(2).YData = r(2, 2, nt);
    end
end