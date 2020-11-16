function energy(planets)
    eng = 0.0
    # kinetic energy
    for p in planets
        eng += 1/2*p.m*norm2(p.v)
    end
    # potential energy
    for j in 1:nplanets
        pj = planets[j]
        for k in j+1:nplanets
            pk = planets[k]
            eng -= G*pj.m*pk.m/sqdist(pj.r, pk.r)
        end
    end
    eng
end
    
function barycenter(m,mTot,x,y,z) # Find Barycenter
    #    m : mass of planet
    # mTot : total mass of system
    # x,y,z: position coordinates of planet
    #
    baryX = (m*x)/mTot
    baryY = (m*y)/mTot
    baryZ = (m*z)/mTot
    return baryX,baryY,baryZ;
end

function momentum(m,x,y,z,vx,vy,vz)
    r = [x,y,z];
    vel = [vx,vy,vz];
    H = m*cross(r,vel);
    return H;
end

