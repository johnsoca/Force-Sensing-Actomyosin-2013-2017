function [NEWxc, NEWyc]=newMotorPos(x1,y1,x2,y2,xc,yc)
% The motor will attach to the filament with plus end at (x1,y1) and minus
% end at (x2,y2)

smin=[(xc-x1)*(x2-x1)+(yc-y1)*(y2-y1)]/[(x2-x1)^2+(y2-y1)^2];

% if smin>0 && smin<1
if smin>1
    smin=1;
elseif smin<0
    smin=0;
end
%     if x2>=x1
        NEWxc=x1+smin*(x2-x1);
%     else
%         NEWxc=x1-smin;
%     end
%     if y2>=y1
        NEWyc=y1+smin*(y2-y1);
%     else
%         NEWyc=y1-smin;
%     end
% else
%     Dp=(xc-x1)^2+(yc-y1)^2;
%     Dm=(xc-x2)^2+(yc-y2)^2;
%     if Dp<Dm %bind to plus end 
%         NEWxc=x1;
%         NEWyc=y1;
%     else
%         NEWxc=x2;
%         NEWyc=y2;
%     end

end