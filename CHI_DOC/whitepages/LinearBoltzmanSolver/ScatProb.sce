clear

N=100;

theta = linspace(0,2*%pi,N);
varphi = linspace(-1*%pi/2,%pi/2,N);

x=[];
y=[];
z=zeros(N*N,N*N);
for i=1:N
    for j=1:N
        x=[x; sin(theta(i))*cos(varphi(N))]
        y=[y; sin(theta(i))*sin(varphi(N))]
        z(i,j) = cos(theta(i));
    end
end

theta0 = linspace(0,%pi,N);
ftheta0 = cos(theta0/2);



scf(0)
clf(0)

surf(x,y,z)
