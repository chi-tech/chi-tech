HarmonicsTest1.luaprint("############### Test of orthogonality")

N = 1000
Nt = N
Nv = 2*N

dtheta = math.pi/Nt
dvaphi = 2.0*math.pi/Nv

ell = 2
m = 1
ellp = 2
mp = 1
intgl = 0.0
for i=0,(Nt-1) do
    theta = dtheta/2.0 + i*dtheta
    for j=0,(Nv-1) do
        vaphi = dvaphi/2.0 + j*dtheta

        Y1 = chiYlm(ell,m,theta,vaphi)
        Y2 = chiYlm(ellp,mp,theta,vaphi)

        intgl = intgl + Y1*Y2*math.sin(theta)*dtheta*dvaphi
    end
end

print("Integral = "..tostring(intgl))
check = 4.0*math.pi/(2*ell+1)
print("Check = "..tostring(check))

print("############### Test of legendre ")
value = 0.0
theta = math.pi/3
vaphi = math.pi/3

x1 = math.sin(theta)*math.cos(vaphi)
y1 = math.sin(theta)*math.sin(vaphi)
z1 = math.cos(theta)

thetap = math.pi/4
vaphip = math.pi/4

x2 = math.sin(thetap)*math.cos(vaphip)
y2 = math.sin(thetap)*math.sin(vaphip)
z2 = math.cos(thetap)

mu = x1*x2 + y1*y2 + z1*z2

for m=(-ell),ell do
    value = value + chiYlm(ell,m,theta,vaphi)*chiYlm(ell,m,thetap,vaphip)
end

print("sum = "..tostring(value))
check = chiLegendre(ell,mu)
print("Check = "..tostring(check))