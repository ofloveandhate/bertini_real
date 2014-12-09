function [] = leap_figure_axes_control(varargin)

global config

% config.mutating_axes = gca;
config.dead_zone = zeros(3,1);

% figure(gcf);

t0 = clock;
while 1
    
    num_pointers = establish_number_pointers();
    
    
    
    
    if num_pointers == 0
        
        if etime(clock, t0)>15
%             display('leap timeout');
            break;
        end
        
        
        continue;
    
    else
%         display('resetting clock');
        t0 = clock;
    end
    
    
    
    if (establish_origin(num_pointers)~=1)
        continue;
    end
    
    
    
    switch num_pointers

        case 1
            orbit_camera();
            
        case 2
            glide();
            
        case 3
            view_pan();
            
        case 4
            
        otherwise
            
            
    end
    
    
    
    
end
    

end




function num_pointers = establish_number_pointers()

% display('establishing number of pointers');
prev_num_pointers = 0;
t_start = clock;
while 1
    
    
    
    
    f = matleap(1);
    
    curr_num_pointers = length(f.pointables);
    if curr_num_pointers==prev_num_pointers%,curr_num_pointers~=0
        if etime(clock,t_start)>1
            break
        end
    else
        t_start = clock;
    end
        
    prev_num_pointers = curr_num_pointers;
    
end

num_pointers = curr_num_pointers;


end



function glide()
global config
% display('gliding');

av_origin = (config.origin(1,:)+ config.origin(2,:))/2;
while 1;
    
%     display('gliding camera');
    current_frame = matleap(1);
    num_fingers = length(current_frame.pointables);
    
    if num_fingers ~= 2
        break
    end
    
    curr_av_pos = (current_frame.pointables(1).position + current_frame.pointables(2).position)/2;
    
    a = (curr_av_pos-av_origin)/100;
%     a = min((current_frame.pointables.position-config.origin),100)/100;
    a(abs(a)<0.10) = 0;
    
    
%     pause
    %get dtheta, dphi from the frame
%     dx = 1e-1;
%     dy = 1e-2;
    dz = 1e-1;
    camdolly(0,0,dz*a(3),'fixtarget');
    drawnow;
end

end


function view_pan()

global config


% av_origin = average_origin();
av_origin = mean(config.origin);%(config.origin(1,:)+ config.origin(2,:)+ config.origin(3,:))/3;
while 1;
    
   
%     display('orbiting camera');
    current_frame = matleap(1);
    num_fingers = length(current_frame.pointables);
    
    if num_fingers == 0
        break
    end
%     current_frame.pointables(1).position
%     current_frame.pointables(2).position
%     cat(1,current_frame.pointables.position)
    
    curr_av_pos = mean(cat(1,current_frame.pointables.position));
    
%     pause
    a = (curr_av_pos-av_origin)/100;
%     a = min((current_frame.pointables.position-config.origin),100)/100;
    a(abs(a)<0.10) = 0;
    

    dtheta = 1;
    dphi = 1;
    campan(dtheta*a(1),dphi*a(2));
    drawnow;
end

end



function orbit_camera()
global config

while 1
    
%     display('orbiting camera');
    current_frame = matleap(1);
    num_fingers = length(current_frame.pointables);
    
    if num_fingers ~= 1
        break
    end
    
    a = (current_frame.pointables.position-config.origin)/100;

	a(abs(a)<0.10) = 0;
    
    
    %get dtheta, dphi from the frame
    dtheta = 10;
    dphi = 5;
    camorbit(dtheta*a(1),dphi*a(2));
    drawnow;
end


end


function success = establish_origin(num_pointers)
global config

% display(sprintf('establishing origin for %i pointers',num_pointers));

config.origin = zeros(num_pointers,3);

t0 = clock;
have_lock = false;

positions = zeros(60,3,num_pointers);
num_consecutive_frames = 0;
num_nonconsecutive_frames = 0;
prevframeid = 0;
while ~have_lock
    
    current_frame = matleap(1);
    
    if current_frame.id == prevframeid
        continue;
    else
        prevframeid= current_frame.id;
    end
    
    num_fingers = length(current_frame.pointables);
    
    if num_fingers ~= num_pointers
        t0 = clock;
        num_consecutive_frames = 1;
        num_nonconsecutive_frames = num_nonconsecutive_frames+1;
%         display('resetting cuz num_pointers');
    else
        num_nonconsecutive_frames = 1;
        num_consecutive_frames = num_consecutive_frames+1;
        
        positions(num_consecutive_frames,:,:) = cat(3,current_frame.pointables.position);
%         positions(1:num_consecutive_frames,:,1)
        
    end
    
    
    
    if and(etime(clock,t0) > 0.5, num_consecutive_frames>=15) %, num_successful_frames
%         positions(:,:,1); % first colon is for all frames.  last
        
        for ii = 1:num_pointers
            config.origin(ii,:) = mean(positions(1:num_consecutive_frames,:,ii));
%             a = mean(positions(1:num_consecutive_frames,:,ii));
%             stats = sqrt(var(positions(1:num_consecutive_frames,:,ii)));
        end
        

        if (1)
            success = true;
            have_lock = true;
%             display('have origin');
        else
%             t0 = t0/2;
%             t0 = 0;
        end
        
       
    end
    
    if num_nonconsecutive_frames > 20
        success = false;
        break
    end
        
end


%pop up a dialog indicating establishing origin.  

%dialog remains up while doing so.


%when done, close dialog and return the origin
end


