function [xhex, yhex, angles]=ellipse(hrad,vrad,x0,y0,numCortex)
   % hrad- horizontal radius of ellipse
   % vrad- vertical radius of ellipse
   % x0,y0- center coordinates
   % xhex- vector of x positions
   % yhex- vector of y positions
   % angles- angles for cortex filaments
   angles=-pi:(2*pi)/numCortex:pi;
   xhex=x0+hrad*cos(angles);
   yhex=y0+vrad*sin(angles);
