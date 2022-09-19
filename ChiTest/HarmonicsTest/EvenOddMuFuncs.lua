Nmu = 5
dtheta = math.pi/2/Nmu
Nazi = 4*3
dvaphi = 2.0*math.pi/Nazi

mus = {}
for k=1,Nmu do
    theta=dtheta/2 + (k-1)*dtheta
    mu = math.cos(theta)
    mus[k] = mu
end

varphis = {}
for k=1,Nazi do
    vaphi = dvaphi/2 + dvaphi*(k-1)
    varphis[k] = vaphi
end

L=5
viz = ""
for ell=0,L do
    for m=-L,L do
        if (math.abs(m)<=ell) then
            odd_in_mu = false
            for kmu,mu in pairs(mus) do
                thetap = math.acos(mu)
                thetam = math.acos(-mu)
                for kv,varphi in pairs(varphis) do
                    Yp = chiYlm(ell,m,thetap,varphi)
                    Ym = chiYlm(ell,m,thetam,varphi)
                    if (math.abs(Yp+Ym)< 1.0e-12 and
                            ((math.abs(Yp)>1.0e-12) and (math.abs(Ym)>1.0e-12))) then
                        odd_in_mu = true
                    end
                    --print(string.format("%3d %3d %7.4f %7.4f %g",ell,m, mu, varphi,Yp))
                    --print(string.format("%3d %3d %7.4f %7.4f %g",ell,m,-mu, varphi,Ym))
                end
            end
            if (odd_in_mu) then
                viz = viz.." 0"
            else
                viz = viz.." x"
            end
        else
            viz = viz.." 0"
        end
    end
    viz = viz.."\n"
end
print("Viz:")
io.write(viz)

--0 0 0 0 0
--0 0 x 0 0
--0 x 0 x x
