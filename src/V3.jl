struct V3{T}
    x::T
    y::T
    z::T
end
@inline Base.zero(::Type{V3{T}}) where T = V3(zero(T), zero(T), zero(T))
@inline Base.zero(::V3{T}) where T = zero(V3{T})
@inline Base.isapprox(a::V3, b::V3; kwargs...) = isapprox(a.x, b.x; kwargs...) && isapprox(a.y, b.y; kwargs...) && isapprox(a.z, b.z; kwargs...)

# +=*/
for OP in [:+, :-]
    @inline @eval function Base.$OP(v::V3, w::V3)
        V3($OP(v.x, w.x), $OP(v.y, w.y), $OP(v.z, w.z))
    end
end

@i @inline @eval function :(+=)(identity)(y!::V3, w::V3)
    y!.x += w.x
    y!.y += w.y
    y!.z += w.z
end

@inline function Base.:(*)(v::V3, x::Real)
    V3(v.x*x, v.y*x, v.z*x)
end
@inline function Base.:(*)(x::Real, v::V3)
    v*x
end

@i @inline function :(+=)(*)(y!::V3, v::V3, x::Real)
    y!.x += v.x * x
    y!.y += v.y * x
    y!.z += v.z * x
end
@i @inline function :(+=)(*)(y!::V3, x::Real, v::V3)
    y! += v*x
end

# norms and distances
@inline function LinearAlgebra.norm2(v::V3)
    sqrt(v.x^2 + v.y^2 + v.z^2)
end
@inline function sqnorm2(v::V3)
    v.x^2 + v.y^2 + v.z^2
end
@inline function sqdistance(r1::V3, r2::V3)
    (r1.x - r2.x)^2 + (r1.y-r2.y)^2 + (r1.z-r2.z)^2
end
@inline function distance(r1::V3, r2::V3)
    sqrt((r1.x - r2.x)^2 + (r1.y-r2.y)^2 + (r1.z-r2.z)^2)
end

@i @inline function :(+=)(sqnorm2)(y!::Real, v::V3)
    y! += v.x^2
    y! += v.y^2
    y! += v.z^2
end

using NiLang.AD
for TP in [:Real, :GVar]
    @eval @i @inline function :(+=)(sqdistance)(y!::$TP, r1::V3{T}, r2::V3{T}) where T
        @routine @invcheckoff begin
            @zeros T dx dy dz
            dx += r1.x - r2.x
            dy += r1.y - r2.y
            dz += r1.z - r2.z
        end
        y! += dx^2
        y! += dy^2
        y! += dz^2
        ~@routine
    end
end