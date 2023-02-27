-- This script tests some basic spherical harmonics functionality


print("Harmonics test script")
quad0_handle = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,32,32)
quad0 = chiGetProductQuadrature(quad0_handle)

quad1_handle = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,4,4)
quad1 = chiGetProductQuadrature(quad1_handle)

--#########################################################
-- This function integrates an angular function F
-- using a quadrature rule
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

--#########################################################
-- This function integrates an angular function F
-- using a simple Riemann-sum with Na azimuthal angles
-- and Np polar angles
function IntegrateFunctionUsingRiemann(F,Na,Np)
    local dvarphi = 2*math.pi/Na
    local dmu = 2.0/Np

    local varphi_min = dvarphi/2
    local mu_min = -1+dmu/2

    local integral = 0.0

    for i=1,Np do
        for j=1,Na do
            local mu = mu_min + (i-1)*dmu
            local theta = math.acos(mu)
            local varphi = dvarphi/2 + (j-1)*dvarphi

            integral = integral + F(theta,varphi)*dmu*dvarphi
        end
    end

    return integral
end

--#########################################################
-- This function integrates an angular function F
-- using a simple Riemann-sum with Na azimuthal angles
-- and Np polar angles with integration limits
function IntegrateFunctionUsingRiemannWL(F,Na,Np, varphi_min, varphi_max, mu_min, mu_max)
    local varphi_range = varphi_max - varphi_min
    local mu_range = mu_max - mu_min
    local dvarphi = varphi_range/Na
    local dmu = mu_range/Np

    local integral = 0.0

    for i=1,Np do
        for j=1,Na do
            local mu = mu_min + dmu/2 + (i-1)*dmu
            local theta = math.acos(mu)
            local varphi = varphi_min + dvarphi/2 + (j-1)*dvarphi

            integral = integral + F(theta,varphi)*dmu*dvarphi
        end
    end

    return integral
end

--#########################################################
--The next three functions are test functions


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

function F3(theta,varphi)
    if (varphi<(7.9*math.pi/8)) then
        return 0.0
    end
    if (varphi>(8.1*math.pi/8)) then
        return 0.0
    end

    return 0.0 + 1.2*math.cos(varphi*40)
end


--#########################################################
-- This function proxy is used to compute Ylm times the function
function_proxy={}
function_proxy["F"] = F0
function_proxy["ell"] = 0
function_proxy["em"] = 0

function F_times_Ylm(theta,varphi)
    local ell = function_proxy.ell
    local em  = function_proxy.em
    return function_proxy.F(theta,varphi)*chiYlm(ell,em,theta,varphi)
end

--#########################################################
-- This function computes the expansion coefficients
-- for a spherical harmonic expansion of the function
-- connected to the function_proxy
function ComputeExpansionCoeffs(expansion_order_L)
    local m=0 --linear moment counter
    local phi_m = {}
    for ell=0,expansion_order_L do
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
            outp = string.format("exp.coeff. m=%3d l=%2d em=%2d %.4e",m,ell,em,value)
            print(outp)
            phi_m[m+1] = value
            m = m + 1
        end--for em
    end--for ell

    return phi_m
end--function



--===============================================
-- This section tests the quadrature integration
print("Testing quadrature integration")
function_proxy.ell = 0
function_proxy.em = 0
integral = IntegrateFunctionUsingQuadrature(F_times_Ylm,quad0)
answer = 4*math.pi
error = 1.0 - integral/answer
print("Integral "..string.format("%.5f",integral)..string.format(" %.5f pi",integral/math.pi))
print("Integral Error "..string.format("%.2e",error))
print("")

--===============================================
-- This section tests the spherical harmonic
-- expansion by computing the expansion
-- coefficients and printing them to temp.txt
print("Testing spherical harmonic expansion")
function_proxy.F = F3
L = 1
phi_m = ComputeExpansionCoeffs(L)
ofile = io.open("temp.txt","w")
for m=1,rawlen(phi_m) do
    if (math.abs(phi_m[m])<1.0e-12) then
        ofile:write("0.0".."\n")
    else
        ofile:write(phi_m[m].."\n")
    end
end
ofile:write("\n")

N=180
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
print("")

--===============================================
-- This section tests the quadrature integration
print("Testing Riemann integration")
function_proxy.ell = 0
function_proxy.em = 0
subdiv = 4
Na = 4*subdiv
Np = 2*subdiv
dvarphi = 2*math.pi/Na
dmu     = 2/Np
max_weight = 0.0
avg_weight = 0.0
counter = 0
for i = 1,Na do
    for j = 1,Np do
        varphi_min = 0.0 + (i-1)*dvarphi
        varphi_max = 0.0 + (i  )*dvarphi

        mu_min = -1.0 + (j-1)*dmu
        mu_max = -1.0 + (j  )*dmu

        integral = IntegrateFunctionUsingRiemannWL(F_times_Ylm,Na*10,Np*10,
                                                   varphi_min,varphi_max,
                                                   mu_min,mu_max)
--         integral = 0.0
        print(string.format("%f %f %f %f %3d %3d %g",varphi_min,varphi_max,
                                                     mu_min,mu_max,
                                                     i,j,integral))
        counter = counter + 1
        max_weight = math.max(max_weight, integral)
        avg_weight = avg_weight + integral
    end
end
avg_weight = avg_weight/counter

print("Max to average: "..string.format("%g", max_weight/avg_weight))


--===============================================
-- This section performs addition theorem checks
print("Addition theorem check")
addition_theorem_works = true
num_angles = rawlen(quad1)
for angle_index=1,num_angles do
    for angle_index_prime=1,num_angles do
        angle_info = quad1[angle_index]
        theta = angle_info.polar
        varphi = angle_info.azimuthal

        omega_x = math.sin(theta)*math.cos(varphi)
        omega_y = math.sin(theta)*math.sin(varphi)
        omega_z = math.cos(theta)

        angle_info = quad1[angle_index_prime]
        theta_prime = angle_info.polar
        varphi_prime = angle_info.azimuthal

        omega_prime_x = math.sin(theta_prime)*math.cos(varphi_prime)
        omega_prime_y = math.sin(theta_prime)*math.sin(varphi_prime)
        omega_prime_z = math.cos(theta_prime)

        mu = omega_x * omega_prime_x +
             omega_y * omega_prime_y +
             omega_z * omega_prime_z;

        for ell=0,1 do
            Pl_mu = chiLegendre(ell,mu)
            sum = 0.0
            for m=-ell,ell do
                val = chiYlm(ell,m,theta,varphi)*
                      chiYlm(ell,m,theta_prime,varphi_prime)
                sum = sum + val
            end

            if (math.abs(sum-Pl_mu) > 1.0e-6) then
                addition_theorem_works = false
            end
        end
    end
end

if (not addition_theorem_works) then
    print("Addition theorem check Failed!")
else
    print("Addition theorem check Passed!")
end
print("")

--===============================================
-- This section performs orthogonality checks
print("Orthogonality checks")
delta_inequal_ell_em_zero = true
delta_equal_ell_em_4pi_2lp1 = true
for ell=0,L do
    for em=-ell,ell do
        for ellp=0,L do
            for emp=-ellp,ellp do

                function YlmYlm(theta,varphi)
                    return chiYlm(ell,em,theta,varphi)*chiYlm(ellp,emp,theta,varphi)
                end

                integral = IntegrateFunctionUsingQuadrature(YlmYlm,quad0)

                if (ell~=ellp or em~=emp) then
                    if (math.abs(integral) > 1.0e-6) then
                        delta_inequal_ell_em_zero = false
                    end
                end
                if (ell==ellp and em==emp) then
                    _4pi_2lp1 = (4.0*math.pi)/(2.0*ell+1.0)
                    if (math.abs(integral-_4pi_2lp1) > 1.0e-6) then
                        delta_equal_ell_em_4pi_2lp1 = false
                        print(integral)
                    end
                end
            end
        end
    end
end

if (not delta_inequal_ell_em_zero) then
    print("Orthogonality check when ell!=ellprime or em!=emprime Failed!")
else
    print("Orthogonality check when ell!=ellprime or em!=emprime Passed!")
end

if (not delta_equal_ell_em_4pi_2lp1) then
    print("Orthogonality check when ell==ellprime and em==emprime Failed!")
else
    print("Orthogonality check when ell==ellprime and em==emprime Passed!")
end
print("")