@i function i_fr!(r::AbstractVector{V3{T}}, v::AbstractVector{V3{T}}, m::AbstractVector; G=G_year_AU, n, dt) where T
    @routine @invcheckoff begin
        nplanets ← length(m)
        θ ← 1/(2-2^(1/3))
        @zeros typeof(dt) halfθdt θdt halfdt_minus_halfθdt dt_minus_twoθdt halfdt
        halfdt += dt / 2
        θdt += dt * θ
        halfθdt += θdt / 2
        halfdt_minus_halfθdt += halfdt - halfθdt
        dt_minus_twoθdt += dt
        dt_minus_twoθdt -= 2 * θdt
    end
    @inbounds @invcheckoff for i=1:n
        a ← zeros(V3{T}, nplanets)
        for j=1:nplanets
            r[j] += halfθdt*v[j]
        end
        @routine i_update_acceleration!(a, r, m, G)
        for j=1:nplanets
            v[j] += θdt*a[j]
        end
        ~@routine
        for j=1:nplanets
            r[j] += halfdt_minus_halfθdt*v[j]
        end
        @routine i_update_acceleration!(a, r, m, G)
        for j=1:nplanets
            v[j] += dt_minus_twoθdt*a[j]
        end
        ~@routine
        for j=1:nplanets
            r[j] += halfdt_minus_halfθdt*v[j]
        end
        @routine i_update_acceleration!(a, r, m, G)
        for j=1:nplanets
            v[j] += θdt*a[j]
        end
        ~@routine
        for j=1:nplanets
            r[j] += halfθdt*v[j]
        end
    end

    ~@routine
end

@i function i_pefrl!(r::AbstractVector{V3{T}}, v::AbstractVector{V3{T}}, m::AbstractVector; G=G_year_AU, n, dt) where T
    @routine @invcheckoff begin
        nplanets ← length(m)
        ξ ← +0.1786178958448091E+00
        λ ← -0.2123418310626054E+00
        χ ← -0.6626458266981849E-01
        @zeros typeof(dt) ξh χh λh halfh_minus_λh h_minus_2χh_minus_2ξh halfh
        halfh += dt / 2
        ξh += ξ * dt
        λh += λ * dt
        χh += χ * dt
        halfh_minus_λh += halfh - λh
        h_minus_2χh_minus_2ξh += dt
        h_minus_2χh_minus_2ξh -= 2 * χh
        h_minus_2χh_minus_2ξh -= 2 * ξh
    end
    @inbounds @invcheckoff for i=1:n
        a ← zeros(V3{T}, nplanets)
        for j=1:nplanets
            r[j] += ξh*v[j]
        end
        @routine i_update_acceleration!(a, r, m, G)
        for j=1:nplanets
            v[j] += halfh_minus_λh*a[j]
        end
        ~@routine
        for j=1:nplanets
            r[j] += χh*v[j]
        end
        @routine i_update_acceleration!(a, r, m, G)
        for j=1:nplanets
            v[j] += λh*a[j]
        end
        ~@routine
        for j=1:nplanets
            r[j] += h_minus_2χh_minus_2ξh*v[j]
        end
        @routine i_update_acceleration!(a, r, m, G)
        for j=1:nplanets
            v[j] += λh*a[j]
        end
        ~@routine
        for j=1:nplanets
            r[j] += χh*v[j]
        end
        @routine i_update_acceleration!(a, r, m, G)
        for j=1:nplanets
            v[j] += halfh_minus_λh*a[j]
        end
        ~@routine
        for j=1:nplanets
            r[j] += ξh*v[j]
        end
    end

    ~@routine
end

