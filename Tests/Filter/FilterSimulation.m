PortLineL = 3e-3;
PortLineW = 1.6e-3;

Whi = 0.2e-3;
Wlo = 3.2e-3;
Lstub = 6.35e-3;
Wstub = 0.238e-3;

L1  = 0.85e-3;
L2  = 3.22e-3;
L3  = 1.54e-3;
L4  = 3.39e-3;
L5  = 2.27e-3;
Lstub1 = 3.45e-3;


a  = traceRectangular('Length',PortLineL,'Width',PortLineW,'Center',[PortLineL/2,0]);
b  = traceRectangular('Length',L1,'Width',Whi,'Center',[PortLineL+L1/2,0]);
c  = traceRectangular('Length',L2,'Width',Wlo,'Center',[PortLineL+L1+L2/2,0]);
l = traceRectangular('Length',Wstub,'Width',Lstub1,'Center',[PortLineL+L1+L2+L3/2,Lstub1/2+Whi/2]);
d  = traceRectangular('Length',L3,'Width',Whi,'Center',[PortLineL+L1+L2+L3/2,0]);
e  = traceRectangular('Length',L4,'Width',Wlo,'Center',[PortLineL+L1+L2+L3+L4/2,0]);
f  = traceRectangular('Length',L5,'Width',Whi,'Center',[PortLineL+L1+L2+L3+L4+L5/2,0]);
g  = traceRectangular('Length',L4,'Width',Wlo,'Center',[PortLineL+L1+L2+L3+L4+L5+L4/2,0]);
m = traceRectangular('Length',Wstub,'Width',Lstub1,'Center',[PortLineL+L1+L2+L3+L4+L5+L4+L3/2,Lstub1/2+Whi/2]);
h  = traceRectangular('Length',L3,'Width',Whi,'Center',[PortLineL+L1+L2+L3+L4+L5+L4+L3/2,0]);
i  = traceRectangular('Length',L2,'Width',Wlo,'Center',[PortLineL+L1+L2+L3+L4+L5+L4+L3+L2/2,0]);
j  = traceRectangular('Length',L1,'Width',Whi,'Center',[PortLineL+L1+L2+L3+L4+L5+L4+L3+L2+L1/2,0]);
k  = traceRectangular('Length',PortLineL,'Width',PortLineW,'Center',[PortLineL+L1+L2+L3+L4+L5+L4+L3+L2+L1+PortLineL/2,0]);
compfiltShape2 = a+b+c+d+e+f+g+h+i+j+k+l+m;
figure;
show(compfiltShape2);

compfilt2 = pcbComponent;
d = dielectric('Teflon');
d.EpsilonR  = 2.2;
d.Thickness = 0.508e-3;
GPL1 = PortLineL+L1+L2+L3+L4+L5+L4+L3+L2+L1+PortLineL;
GPW  = 20e-3;
gnd  = traceRectangular('Length',GPL1,'Width',GPW/2,'Center',[GPL1/2,0]);
compfilt2.BoardThickness = 0.508e-3;
compfilt2.Layers         = {compfiltShape2,d,gnd};
compfilt2.BoardShape     = gnd;
compfilt2.FeedDiameter   = PortLineW/2;
compfilt2.FeedLocations  = [0,0,1,3;GPL1,0,1,3];
compfilt2.ViaLocations   = [PortLineL+L1+L2+L3/2,Lstub1+Whi/2-0.2e-3,1,3;PortLineL+L1+L2+L3+L4+L5+L4+L3/2,Lstub1+Whi/2-0.2e-3,1,3];
compfilt2.ViaDiameter    = Wstub/2;
figure;
show(compfilt2);

spar5 = sparameters(compfilt2,linspace(0.5e9,14e9,41));
figure;
rfplot(spar5);