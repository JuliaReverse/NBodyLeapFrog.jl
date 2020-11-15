"""
Created on Mon Mar 21 20:58:03 2016
Computational Physics
@author: joshualamstein

This script simulates the movement of the planets and the Sun in our solar system. Using a 2nd order leap frog method, 
the simulation conserves energy. To keep track of any mistakes or uncertainties, the change in energy and change in momentum
of the solar system are plotted. Small uncertainties, on the order of 10^-7, are expected from roundoff error. 

If you change the variable checkMaximallyPacked=0 to checkMaximallyPacked=1, 
the simulations adds a planetX between Mars and Jupiter, showing that it is possible
to add a planet with mass similar to Earth to our solar system and still maintain stable orbit. 
Between Mars and Jupiter is ideal since that's where the asteroid belt is. 

---------------------------------------
Translated to Julia by: GiggleLiu, Nov 15, 2020
"""

# ==============================================================================
#  SET UP
# ==============================================================================

solar = 1.988544*10^30 # in kg
year = 3.154*10^7 #year in seconds
jup = 1898.13*10^24 #in kg
earth = 5.97219e24 #in kg
AU = 1.496*10^11 #in m

dayToYear = 365.25; 

# Convert gravitational constant to astronomical units (AU) and solar masses

G_old = 6.67259e-11 # in m^3/(kg-s^2)
newG = G_old*(1/AU)^3/(1/solar*(1/(year))^2)
G = newG # Gravitational Constant

N = 4000; #Iterations
dt = .01; #Time step (years)
checkMaximallyPacked = 0; # CHANGE TO 1 TO ADD A PLANET BETWEEN MARS AND JUPITER AND MAINTAIN STABLE ORBITS

if (checkMaximallyPacked ==1)
    planets =11; #Includes Pluto, remembering the good ol' days. It also includes the Sun.
    x,y, z,vx ,vy ,vz = zeros((planets,N+1)),zeros((planets,N+1)),zeros((planets,N+1)),
        zeros((planets,N+1)),zeros((planets,N+1)),zeros((planets,N+1));
    vxPacked = vxPacked*dayToYear;
    vyPacked = vyPacked*dayToYear;
    vzPacked = vzPacked*dayToYear 
    m = mPacked;
    for i in 1:planets
        #These arrays come the starData.py


        x[i,1] = xPacked[i][1];
        y[i,1] = yPacked[i][1];
        z[i,1] = zPacked[i][1];

        vx[i,1] = vxPacked[i][1];
        vy[i,1] = vyPacked[i][1];
        vz[i,1] = vzPacked[i][1];
    end
else
    planets =10; #Includes Pluto, remembering the good ol' days. It also includes the Sun. 
    x,y, z,vx ,vy ,vz = zeros((planets,N+1)),zeros((planets,N+1)),zeros((planets,N+1)),
    zeros((planets,N+1)),zeros((planets,N+1)),zeros((planets,N+1));
    vxSet = vxSet*dayToYear;
    vySet = vySet*dayToYear;
    vzSet = vzSet*dayToYear
    for i in 1:planets
        #These arrays come the starData.py

        x[i,1] = xSet[i][1];
        y[i,1] = ySet[i][1];
        z[i,1] = zSet[i][1];
        vx[i,1] = vxSet[i][1];
        vy[i,1] = vySet[i][1];
        vz[i,1] = vzSet[i][1];  
    end
end


timeSpace = 0:dt:dt*N
timePeriod = 0:dt:dt*10
ChangeInEnergy = zeros(N); # Record change in energy to check validity of code
h = zeros((N,planets,3)) # variable for momentum
ChangeInH = zeros(N); # change in momentum
totalH = zeros((N,3)) # total momentum
hChange = zeros(N) #placeholder for momentum to check code

energy = zeros(N); 

mag = zeros((planets,planets)) # distances between planets
#acceleration
axNew,ayNew,azNew = zeros((planets,planets)),zeros((planets,planets)),zeros((planets,planets))
ax = zeros((planets,N))
ay = zeros((planets,N))
az = zeros((planets,N))
axPull, ayPull,azPull = zeros((planets,1)),zeros((planets,1)),zeros((planets,1))

#position and velocity at Cartesian vector components

vBary= zeros((N,3)) # Velocity barycenter
rBary = zeros((N,3)) #barycenter position

  



# ==============================================================================
#  function definitions
# ==============================================================================
    
function magnitude(x,y,z)
    # x,y,z : coordinates of planet
    return sqrt(x^2 + y^2 +z^2)
end
    
function momentum(m,x,y,z,vx,vy,vz)
    r = [x,y,z];
    vel = [vx,vy,vz];
    H = m*cross(r,vel);
    return H;
end

# Update position of planet with Leap frog method, 
# 2*dt is the length of the time step for leap frog
function position(xOld,yOld,zOld,vx,vy,vz,dt)
    xNew = xOld+2*dt*vx
    yNew = yOld+2*dt*vy
    zNew = zOld+2*dt*vz
    return xNew,yNew,zNew
end

function positionEuler(x,y,z,vx,vy,vz,dt) # Get position using the Euler method
    xNew = x + vx*dt;
    yNew = y + vy*dt;
    zNew = z + vz*dt;
    return xNew,yNew,zNew;
end
    
function acceleration(x,y,z,xPull,yPull,zPull, mag, G,mPull) # Get acceleration of celestial body
    ax = G * mPull * (xPull-x)/mag^3;
    ay = G * mPull * (yPull-y)/mag^3;
    az = G * mPull * (zPull-z)/mag^3;
    return ax,ay,az
end
    
function velocity(vxOld,vyOld,vzOld,ax,ay,az,dt) # Get velocity of planet with leap frog method
    vxNew = vxOld + 2*dt*ax
    vyNew = vyOld + 2*dt*ay
    vzNew = vzOld + 2*dt*az

    return vxNew,vyNew,vzNew;
end
    
function velocityEuler(vx,vy,vz,ax,ay,az,dt) #Euler method
    vxNew = vx + ax*dt;
    vyNew = vy + ay*dt;
    vzNew = vz + az*dt;
    
    return vxNew,vyNew,vzNew;
end
    
function kinetic(m,vx,vy,vz) # Kinetic Energy
    KE = 1/2 *m* (vx^2+vy^2+vz^2);
    return KE;
end
    
function potential(G,m,mPull,mag) # Gravitational Potential
    #     G : Gravitational constant
    #     m : mass of planet
    # mPull : mass of other object (probably the Sun)
    #   mag : magnitude displacement
    #
    U = -G*m*mPull/mag;
    return U;
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

function nOrbit(x,y,z,xPull,yPull,zPull,vx,vy,vz,vxPull,vyPull,vzPull,G,m,mPull,n,planets,dt) # this is the bulk of the code
    #                  x,y,z : inital position coordinates of astronomical body
    #    xPull, yPull, zPull : initial position coordinates of other astronomical body interacting with astro body #1
    #               vx,vy,vz : initial velocity components of astro body #1
    # vxPull, vyPull, vzPull : initial velocity components of astro body #1
    #                      G : Gravitational constant
    #                      m : mass of planet
    #                  mPull : mass of other object (probably the Sun)
    #                      n : Number of iterations
    #                     dt : Time step (years)
    mTot = 0;
    rAdd = 0;
    vAdd = 0;
    for j in 1:planets
         mTot += m[j] # total mass of solar system
    end
                                                      
    for i in 1:n
        
        
        # ==============================================================================
        # acceleration
        # ==============================================================================
            
        #Get distance between all included celestial bodies
        for j in 1:planets
            for k in j:planets
                mag[j,k] = magnitude(x[j,i]-xPull[k,i],y[j,i]-yPull[k,i],z[j,i]-zPull[k,i]);
            end
            if j > 0
                for p in 1:j
                    mag[j,p] = mag[p,j]; #symmetry, distance between one planet is the same as the distance from the other
                end
            end
        end

        # Find acceleration due to gravity
        for j in 1:planets
            for k in 1:planets
                if j != k
                     axNew[j,k],ayNew[j,k],azNew[j,k] = acceleration(x[j,i],y[j,i],z[j,i],
                     xPull[k,i], yPull[k,i],zPull[k,i],mag[j,k],G,mPull[k]);
                end
            end
            sumAxNew,sumAyNew,sumAzNew = sum(axNew[j]),sum(ayNew[j]),sum(azNew[j])                                       
            ax[j,i],ay[j,i],az[j,i] = sumAxNew,sumAyNew,sumAzNew
        end

        # ==============================================================================
        #   velocity  
        # ==============================================================================

        #At the beginning of the array, you can't use leap frog, so it suffices to use the 1st order Euler method
        if i==1
            for j in 1:planets
                vxNew,vyNew,vzNew = velocityEuler(vx[j,i],vy[j,i],vz[j,i],ax[j,i],ay[j,i],az[j,i],dt)
                vx[j,i+1],vy[j,i+1],vz[j,i+1] = vxNew,vyNew,vzNew;
            end
        else # Use leap frog method to update velocity
            for j in 1:planets
                vxNew,vyNew,vzNew = velocity(vx[j,i-1],vy[j,i-1],vz[j,i-1],ax[j,i],ay[j,i],az[j,i],dt);
                vx[j,i+1],vy[j,i+1],vz[j,i+1] = vxNew,vyNew,vzNew;
            end
        end
                 
            
        # ==============================================================================
        # Energy        
        # ==============================================================================
        # Because energy is conserved,
        # the change in energy should be due to round off error in the simulation. 

        KEhold = 0; #Kinetic energy placeholder
        PEhold = 0; #Potential energy placeholder

        for j in 1:planets
            KEhold += kinetic(m[j],vx[j,i],vy[j,i],vz[j,i])
        end
        for j in 1:planets
            for k in j:planets
                if j != k
                    PEhold += potential(G,m[j],mPull[k], mag[j,k])
                end
            end
        end
                    
        energy[i] = KEhold + PEhold;

        ChangeInEnergy[i] = abs(energy[i])-abs(energy[1]);


        ## ==============================================================================
        ## check barycenter which should be 0
        ## ==============================================================================
        for j in 1:planets
            rAdd = barycenter(m[j],mTot,x[j,i],y[j,i],z[j,i],)
            rBary[i,:] .= rBary[i] .+ rAdd
            vAdd = barycenter(m[j],mTot,vx[j,i],vy[j,i],vz[j,i],)
            vBary[i,:] .= vBary[i] .+ vAdd        
        end
        rAdd = 0;
        vAdd= 0;

        # ==============================================================================
        #  new position
        # ==============================================================================
            
        # Get position with Euler method for i==0      
        if i ==1
            for j in 1:planets
                xNew,yNew,zNew = positionEuler(x[j,i],y[j,i],z[j,i],vx[j,i],
                                      vy[j,i],vz[j,i],dt) ;
                x[j,i+1],y[j,i+1],z[j,i+1] = xNew,yNew,zNew;
                xPull[j,i+1],yPull[j,i+1],zPull[j,i+1] = xNew,yNew,zNew;
            end
            # Update position with leap frog method
        else
            for j in 1:planets
                xNew,yNew,zNew = position(x[j,i-1],y[j,i-1],z[j,i-1],vx[j,i],
                                      vy[j,i],vz[j,i],dt);
                x[j,i+1],y[j,i+1],z[j,i+1] = xNew,yNew,zNew;
                xPull[j,i+1],yPull[j,i+1],zPull[j,i+1] = xNew,yNew,zNew;
            end
        end

        ## ==============================================================================
        ## momentum (no mass)
        ## ==============================================================================
        for j in 1:planets
            h[i,j,:] .= momentum(m[j],x[j,i],y[j,i],z[j,i],vx[j,i],vy[j,i],vz[j,i])
            totalH[i,:] .= dropdims(sum(h[i,:,:],dims = 1), dims=1)
        end
    end
    hTotalMag = norm(totalH)
    hTotComb = [norm(totalH[i,:]) for i=1:size(totalH, 1)]
    ChangeInH = hTotComb .- hTotComb[1]
    hRatio = ChangeInH ./ hTotalMag
    
    # ChangeInEnergy,hMag,hPullMag,bary,vbary

    return x,y,z,vx,vy,vz,ax,ay,az,energy, ChangeInEnergy, totalH,ChangeInH,hRatio, hTotalMag,rBary,vBary;
end


# ==============================================================================
#  RUN SCRIPT
# ==============================================================================
    
a,b,c,va,vb,vc,aa,ab,ac,energy, ChangeInE,totalH, ChangeInH,hRatio,hTotalMag, rBary,vBary = nOrbit(x,y,z,x,y,z,vx,vy,vz,vx,vy,vz,G,m,m,N,planets,dt);

plt = plot([], [])
for i in 1:planets
    plot!(plt,a[i,:],b[i,:])
end
