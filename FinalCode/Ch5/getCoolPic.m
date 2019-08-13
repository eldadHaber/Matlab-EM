
absb = Af*(abs(diff(b,[],2)));

% Mv = moviein( size(absb,2));
% get video writer object
writeObj = VideoWriter('currentDensity.avi');
writeObj.FrameRate = 10;
open(writeObj);


for i=1:size(absb,2),
    bb = reshape(absb(:,i),n1,n2,n3); 
    %trs = 1e-6; %max(bb(:))/500
    %clf;volview(flipdim(bb,3),'isovalue',trs); 
    slice(interpn(bb,2),65,65,1);
    shading flat
    axis off tight
    view(48,32)
    figure(1)
    frame = getframe;
    writeVideo(writeObj,frame);
    
    
    %pause; 
end

close(writeObj);