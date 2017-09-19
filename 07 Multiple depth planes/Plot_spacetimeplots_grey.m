

depths = [25:-5:0,-10,-20];
N_depths = length(depths);
ground = cell(N_depths,1);

for nd = 1:N_depths
    
    screen_distance = 20.0; % cm
    screen_height = 26.3; % cm
    
    grating_wavelength = 10; % cm
    N_periods = 40; % That are drawn on the floor
    
    grating_tf = 2; % Temporal frequency of the grating
    grating_speed = grating_tf*grating_wavelength; % cm/s
    
    frame_rate = 200; % Hz
    dx = grating_speed / frame_rate; % grating movement per frame, cm/frame
    
    stimulus_duration = 4.0; % s
    N_frames =  stimulus_duration * frame_rate;
    
    step_height = depths(nd); % cm
    step_width = 1 * grating_wavelength;
    
    step_start = grating_speed*stimulus_duration; % to end arriving at the bootom of screen
    step_start = step_start - step_width - grating_speed*0.5; % to end with grating passing out of the screen 0.5s before the end when step height is zero
    step_start = step_start - step_height*screen_distance/screen_height; % Correct for height of step
    
    ground{nd} = zeros(480, N_frames)+0;
    ground{nd}(240:480,:) = 0.5;
    
    for nf = 1:N_frames
        % Draw ground
        x_start = mod(-dx*(nf-1), grating_wavelength);
        x = x_start-grating_wavelength:grating_wavelength/2:N_periods*grating_wavelength; % horizontal distance behind screen
        y = screen_height-(screen_distance*screen_height)./(screen_distance + x); % cm
        y = round(y *480/screen_height); % pixels
        
        for n = 1:2:N_periods-2
            id = y(n)+1:y(n+2)+1;
            pos_id = find(id > 0);
            ground{nd}(id(pos_id),nf) = 0.5;
        end
        
        % Draw step
        x2_start = step_start            - (nf-1)*dx;
        x2_end   = step_start+step_width - (nf-1)*dx;
        x2 = x2_start:grating_wavelength/2:x2_end+0.0001; % the 0.0001 helps correct for rounding errors
        for n = 1:length(x2)
            if (x2(n) < -screen_distance+1)
                x2(n) = -screen_distance+1;
            end
        end
        y2 = screen_height-(screen_distance*(screen_height-step_height))./(screen_distance + x2); % cm
        y2 = round(y2 *480/screen_height); % pixels
        for n = 1:length(y2)
            if (y2(n) < 0)
                y2(n) = 0;
            end
        end
        for n = 1:2:length(y2)-2
            ground{nd}(y2(n)+1:y2(n+1)+1,nf) = 0; % 0.4
            ground{nd}(y2(n+1)+1:y2(n+2)+1,nf) = 1; % 0.9
        end
    end
    
    figure(100)
    subplot(2,4,nd)
    hold off
    mesh(ground{nd})
    axis([0 N_frames 0 480 0 2])
    view(2)
    colormap('bone')
    set(gca,'xtick',0:200:N_frames,'xticklabel',{'0' '1' '2' '3' '4'});
        
   
end

for nl = 1:N_depths
   subplot(2,4,nl)
   xlabel('Time (s)')
   ylabel('Y-axis pixels')
end

save spacetimeplots_grey.mat ground





