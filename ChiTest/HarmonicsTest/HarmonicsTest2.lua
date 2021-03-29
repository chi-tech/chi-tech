print("Hello")
quad0_handle = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,32,32)
quad0 = chiGetProductQuadrature(quad0_handle)

--===============================================
function IntegrateFunctionUsingQuadrature(F,quad)
    local integral = 0.0

    local num_angles = rawlen(quad0)
    for angle_index=1,num_angles do
        angle_info = quad0[angle_index]
        theta = angle_info.polar
        varphi = angle_info.azimuthal
        integral = integral + angle_info.weight*F(theta,varphi)
    end

    return integral
end

function IntegrateFunctionUsingRiemann(F,Na,Np)
    dvarphi = 2*math.pi/Na
    dmu = 2.0/Np

    varphi_min = dvarphi/2
    mu_min = -1+dmu/2

    integral = 0.0

    for i=1,Np do
        for j=1,Na do
            mu = mu_min + (i-1)*dmu
            theta = math.acos(mu)
            varphi = dvarphi/2 + (j-1)*dvarphi

            integral = integral + F(theta,varphi)*dmu*dvarphi
        end
    end

    return integral
end

function F0(theta,varphi)
    return 1.0
end

function F1(theta,varphi)
    return 1.0 + 0.1*math.cos(varphi*4)
end

function F2(theta,varphi)
    if (varphi<(7*math.pi/8)) then
        return 0.2
    end
    if (varphi>(9*math.pi/8)) then
        return 0.2
    end

    return 0.2 + 1.2*math.cos(varphi*4)
end

function_proxy={}
function_proxy["F"] = F2
function_proxy["ell"] = 0
function_proxy["em"] = 0

function F_times_Ylm(theta,varphi)
    local ell = function_proxy.ell
    local em  = function_proxy.em
    return function_proxy.F(theta,varphi)*chiYlm(ell,em,theta,varphi)
end

L = 7
m=0
phi_m = {}
for ell=0,L do
    for em=-ell,ell do
        function_proxy.ell = ell
        function_proxy.em = em

        function T(theta,varphi)
            local ell = function_proxy.ell
            local em  = function_proxy.em
            return function_proxy.F(theta,varphi)*chiYlm(ell,em,theta,varphi)
        end

        value = IntegrateFunctionUsingQuadrature(T,quad0)
--         value = IntegrateFunctionUsingRiemann(T,100,100)
        outp = string.format("m=%3d l=%2d m=%2d %.4e",m,ell,em,value)
        print(outp)
        phi_m[m+1] = value
        m = m + 1
    end
end

--===============================================
integral = IntegrateFunctionUsingQuadrature(F_times_Ylm,quad0)
answer = 4*math.pi
error = 1.0 - integral/answer
print("Integral "..string.format("%.5f",integral)..string.format("%.5f pi",integral/math.pi))
print("Integral Error "..string.format("%.2e",error))

ofile = io.open("temp.txt","w")
for m=1,rawlen(phi_m) do
    if (math.abs(phi_m[m])<1.0e-12) then
        ofile:write("0.0".."\n")
    else
        ofile:write(phi_m[m].."\n")
    end
end
ofile:write("\n")

N=200
dvarphi = 2*math.pi/N
dtheta = math.pi/N

for i=1,N do
    varphi= dvarphi/2 + (i-1)*dvarphi

    func_val = function_proxy.F(0,varphi)

    approx = 0.0
    m=0
    for ell=0,L do
        for em=-ell,ell do
            mu = 0.0
            theta = math.acos(mu)
            approx = approx + ((2.0*ell+1.0)/4.0/math.pi)*
                              phi_m[m+1]*chiYlm(ell,em,theta,varphi)
            m = m + 1
        end
    end

    outp = string.format("%.4f ",varphi)
    outp = outp..string.format("%.4e ",func_val)
    outp = outp..string.format("%.10e",approx)
    outp = outp.."\n"

    ofile:write(outp)
end

ofile:close(ofile)