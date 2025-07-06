function animation(h,index,filename)
    
    % Inputs:
    %   h        - Handle to the figure to capture
    %   index    - Frame index
    %   filename - Name of the output GIF file

    % Output:
    %   GIF file.

    % Author: Jonatan Ismael Eisermann  
    % Date: July 6, 2025. 

    drawnow
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    if index == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end

end