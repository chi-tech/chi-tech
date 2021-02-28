quad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,18,90)

quad_raw = chiGetProductQuadrature(quad)

N = rawlen(quad_raw)
print("Number of angles = "..tostring(N)  )

total_weight = 0.0
for n=1,N do
    print(string.format("Weight=%9.5f Polar=%9.5f Azimuthal=%9.5f",
          quad_raw[n].weight,
          quad_raw[n].polar,
          quad_raw[n].azimuthal))
    total_weight = total_weight + quad_raw[n].weight
end

print("Total weight = "..tostring(total_weight))

--===================================== Compute angular expansion
L = 7
-- function F(theta,varphi)
--     return 1.0+0.1*math.cos(varphi*4.0)
-- end

function F(theta,varphi)
    if (varphi<(7*math.pi/8)) then
        return 0.0
    end
    if (varphi>(9*math.pi/8)) then
            return 0.0
    end

    return 1.2*math.cos(varphi*4)+0.0
end

--===================================== Compute expansion coeffs
f_mstar = {}
C_mstar = {}
mstar=0
for ell=0,L do
    for m=(-ell),ell do
        mstar = mstar + 1
        f_mstar[mstar] = 0.0
        C_mstar[mstar] = (2*ell+1.0)/4.0/math.pi

        for n=1,N do
            w = quad_raw[n].weight
            theta = quad_raw[n].polar
            phi = quad_raw[n].azimuthal
            Ylm = chiYlm(ell,m,theta,phi)

            if ((phi>0.5*math.pi) and ((phi<0.5*3*math.pi))) then
                f_mstar[mstar] = f_mstar[mstar] + w*F(theta,phi)*Ylm
            end
        end
    end
end

--===================================== Compute a sample angular distro
Na = 360
dphi = 2.0*math.pi/Na
for a=0,Na do
    theta = 0.5*math.pi
    phi = dphi*a

    expansion = 0.0

    mstar=0
    for ell=0,L do
        for m=(-ell),ell do
            mstar = mstar + 1
            Ylm = chiYlm(ell,m,theta,phi)
            --print(ell,m,C_mstar[mstar],f_mstar[mstar],Ylm)
            expansion = expansion + C_mstar[mstar]*f_mstar[mstar]*Ylm
        end
    end

    print(string.format("%9.5f %9.5f %9.5f",phi,expansion,F(theta,phi)))
end

