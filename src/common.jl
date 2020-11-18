using .Bodies: G_year_AU, Body

struct LeapFrog end
struct PEFRL end
struct FR end

function leapfrog(::Type{LF}, planets::AbstractVector{Body{T}}; G=G_year_AU, n, dt, keep_history=false) where {T, LF}
    nplanets = length(planets)
    m = getfield.(planets, :m)
    @inbounds if keep_history
        v = zeros(V3{T}, (nplanets,n+1))
        r = zeros(V3{T}, (nplanets,n+1))
        v[:,1] .= getfield.(planets, :v)
        r[:,1] .= getfield.(planets, :r)
    else
        v = zeros(V3{T}, nplanets)
        r = zeros(V3{T}, nplanets)
        v .= getfield.(planets, :v)
        r .= getfield.(planets, :r)
    end
    leapfrog_looper!(LF, r, v, m; G=G, dt=dt, n=n)
end

function leapfrog_looper!(::Type{LF}, r::AbstractVector{V3{T}}, v::AbstractVector{V3{T}},
         m::AbstractVector; G=G_year_AU, n, dt) where {T, LF}
    for i=1:n
        steper(LF)(r, v, m; G=G, dt=dt)
    end
    return r, v
end

function leapfrog_looper!(::Type{LF}, r::AbstractMatrix{V3{T}}, v::AbstractMatrix{V3{T}},
         m::AbstractVector; G=G_year_AU, n, dt) where {T,LF}
    for i=1:n
        v[:,i+1] .= view(v,:,i)
        r[:,i+1] .= view(r,:,i)
        steper(LF)(view(r,:,i+1), view(v,:,i+1), m; G=G, dt=dt)
    end
    return r, v
end

@inline function acceleration(ra, rb, mb, G)
    d = distance(ra, rb)
    (G*mb/d^3) * (rb - ra)
end

function update_acceleration!(a::AbstractVector{V3{T}}, r, m, G) where T
    @assert length(a) == length(m)
    @assert length(a) == length(r)
    @inbounds for j=1:length(m)
        a[j] = zero(V3{T})
        for k=1:length(m)
            if j != k
                a[j] += acceleration(r[j], r[k], m[k], G)
            end
        end
    end
    return a
end

