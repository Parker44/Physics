function [t, r] = toomre(tmax, level, cores_m, cores_r0, cores_v0, cores_ns)
    
    % Discretize continuum domain
    nt = 2^(level) + 1;
    deltat = tmax / (nt - 1);
    t = zeros(1, nt);
    for i=2:nt
       t(i) = t(i-1) + deltat;
    end
    
    n_cores = length(cores_m);
    plt = gobjects(1, n_cores);
   
    % Initialize 
    r = zeros(n_cores, 1 + cores_ns, 3, nt);
    v_0 = zeros(n_cores, 1 + cores_ns, 3);
    
    
    for i=1:n_cores
        r(i, 1, :, 1) = cores_r0(i, :);
        v_0(i, 1, :) = cores_v0(i, :);
        
        r(i, 1, :, 2) = r(i, 1, :, 1) + deltat*v_0(i, 1, :);
        
        for j=2:cores_ns+1
            r_min = 0.03;  % minimum star orbital radius
            r_max = 0.2; % maximum star orbital radius
            % Set random radius and theta
            d = (r_max-r_min)*rand(1,1) + r_min;
            theta = (2*pi)*rand(1,1);
            r_orbit = [d*cos(theta) d*sin(theta) 0];
            
            r(i, j, :, 1) = reshape(r(i, 1, :, 1), [1,3]) + r_orbit;
            v_0(i, j, :) = reshape(v_0(i, 1, :), [1,3]) + sqrt(cores_m(i))*norm(r_orbit)^(-1/2)*[-sin(theta) cos(theta) 0];
            r(i, j, :, 2) = reshape(r(i, j, :, 1), [1,3]) + deltat*reshape(v_0(i, j, :), [1,3]);
        end
    end
    
    a = nbodyaccn(cores_m, n_cores, cores_ns, r(:, :, :, 1));
    
    for i=1:n_cores
        r(i, 1, :, 2) = r(i, 1, :, 1) + deltat*v_0(i, 1, :) + 0.5 * deltat^2 * a(i, 1, :);
        
        for j=2:cores_ns+1
            r(i, j, :, 2) = r(i, j, :, 1) + deltat*v_0(i, j, :) + 0.5 * deltat^2 * a(i, j, :);
        end
    end
    %-----------------------------------------------------------
    % Set plotenable to non-zero/zero to enable/disable plotting.
    %-----------------------------------------------------------
    plotenable = 1;
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

    % Particles have a (marker) size of 2
    bodysize = 2;
    % they are yellow
    bodycolors = ['y', 'g', 'm', 'c'];
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
    avifilename = 'toomre.avi';

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

        % Setup
        axis square;
        xlim([-1, 1]);
        ylim([-1, 1]);
        set(gca, 'Color', bgcolor);
        set(gcf, 'Visible', 'off');
        set(gcf,'Position',[100 100 250 250])
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        
        for i=1:n_cores
            plt(i) = plot(r(i, :, 1, 1), r(i, :, 2, 1), 'Marker', bodymarker, 'MarkerSize', bodysize, ...
                'MarkerEdgeColor', bodycolors(i), 'MarkerFaceColor', bodycolors(i), 'LineStyle', 'none');
        end
    end
    
    % plot
    for n=2:nt-1
        
        if plotenable
            
            for i=1:n_cores
                % Draw the particles
                plt(i).XData = r(i, :, 1, n);
                plt(i).YData = r(i, :, 2, n);
            end
                
            % Force update of figure window.
            % drawnow;

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
        
        % Calculate acceleration matrix (all particles)
        a = nbodyaccn(cores_m, n_cores, cores_ns, r(:, :, :, n));
        for i=1:n_cores
            for j=1:cores_ns+1
                % Calculate next position
                r(i, j, :, n+1) = reshape(2*r(i, j, :, n), [1,3]) - reshape(r(i, j, :, n-1), [1,3]) + (deltat^2)*reshape(a(i,j, :), [1,3]);
            end 
        end
    end
    %Plot last point
    if plotenable
        plt(i).XData = r(i, j, 1, nt);
        plt(i).YData = r(i, j, 2, nt);
    end
end