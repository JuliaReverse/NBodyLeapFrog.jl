using .Bodies: V3, distance, mass_solar
#  unit of time -> year, unit of space -> AU
const year = 3.154e7 #year in seconds
const AU = 1.496e11 #in m

const G_standard = 6.67259e-11 # in m^3/(kg-s^2)
const G_year_AU = G_standard*(1/AU)^3/(1/mass_solar*(1/year)^2)

function acceleration(ra, rb, mb, G) # Get acceleration of celestial body
    d = distance(ra, rb)
    (G*mb/d^3) * (rb - ra)
end

function nOrbit(planets; G=G_year_AU, n, dt) # this is the bulk of the code
    nplanets = length(planets)
    a = zeros(V3{Float64}, (nplanets,n))
    v = zeros(V3{Float64}, (nplanets,n+1))
    r = zeros(V3{Float64}, (nplanets,n+1))
    v[:,1] .= getfield.(planets, :v)
    r[:,1] .= getfield.(planets, :r)
    for i in 1:n
        # compute acceleration
        for j=1:nplanets
            for k=1:nplanets
                if j != k
                    a[j,i] += acceleration(r[j,i], r[k,i], planets[k].m, G)
                end
            end
        end

        # update velocity  
        if i==1
            for j in 1:nplanets
                v[j,i+1] = dt*a[j,i] + v[j,i]
            end
        else # Use leap frog method to update velocity
            for j in 1:nplanets
                v[j,i+1] = 2dt*a[j,i] + v[j,i-1]
            end
        end
                 
        # update position
        if i ==1
            # Get position with Euler method for i==0      
            for j in 1:nplanets
                r[j,i+1] = dt*v[j,i] + r[j,i]
            end
        else
            # Update position with leap frog method
            for j in 1:nplanets
                r[j,i+1] = 2dt*v[j,i] + r[j,i-1]
            end
        end
    end
    return r, v, a
end