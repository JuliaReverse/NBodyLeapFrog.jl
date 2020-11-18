
function leapfrog_step!(r::AbstractVecOrMat{V3{T}}, v::AbstractVecOrMat{V3{T}},
         m::AbstractVector; G, dt) where T
    @assert length(r) == length(m) == length(v)
    a = zeros(V3{T}, length(m))
    @inbounds for j=1:length(m)
        r[j] = dt/2*v[j] + r[j]
    end
    update_acceleration!(a, r, m, G)
    @inbounds for j=1:length(m)
        v[j] = dt*a[j] + v[j]
        r[j] = dt/2*v[j] + r[j]
    end
    return r, v
end

function leapfrog!(r::AbstractVecOrMat{V3{T}}, v::AbstractVecOrMat{V3{T}},
         m::AbstractVector; G=G_year_AU, n, dt) where T
    leapfrog_looper!(LeapFrog, r, v, m; G=G, n=n, dt=dt)
end

function leapfrog(planets::AbstractVector{Body{T}}; G=G_year_AU, n, dt, keep_history=false) where {T}
    leapfrog(LeapFrog, planets; G=G, n=n, dt=dt, keep_history=keep_history)
end

steper(::Type{LeapFrog}) = leapfrog_step!