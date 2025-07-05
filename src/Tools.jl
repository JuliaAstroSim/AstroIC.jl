"""
    function addpos(data::Array, pos::AbstractPoint)
    function addpos(data::StructArray, pos::AbstractPoint)

Add `pos` to `:Pos` of all particles
"""
function addpos(data::Array, pos::AbstractPoint)
    for i in 1:length(data)
        data[i] = setproperties!!(data[i], Pos = data[i].Pos + pos)
    end
end

function addpos(data::StructArray, pos::AbstractPoint)
    Pos = data.Pos
    for i in 1:length(data)
        Pos[i] += pos
    end
end

"""
    function addvel(data::Array, vel::AbstractPoint)
    function addvel(data::StructArray, vel::AbstractPoint)

Add `vel` to `:Vel` of all particles
"""
function addvel(data::Array, vel::AbstractPoint)
    for i in 1:length(data)
        data[i] = setproperties!!(data[i], Vel = data[i].Vel + vel)
    end
end

function addvel(data::StructArray, vel::AbstractPoint)
    Vel = data.Vel
    for i in 1:length(data)
        Vel[i] += vel
    end
end

"""
    function setpos(data::Union{Array, StructArray}, pos::AbstractPoint)

Set system center (middle value) to `pos`
"""
function setpos(data::Union{Array, StructArray}, pos::AbstractPoint)
    p0 = median(data, :Pos)
    addpos(data, pos - p0)
end

"""
    function setvel(data::Union{Array, StructArray}, vel::AbstractPoint)

Set system velocity (mass weighted average) to `vel`
"""
function setvel(data::Union{Array, StructArray}, vel::AbstractPoint)
    v0 = averagebymass(data, :Vel)
    addvel(data, vel - v0)
end


"""
$(TYPEDSIGNATURES)
Generate `PVector2D` from radius
"""
function rand_pos_2d(r::Number)
    gamma = rand()

    x = r * cos(2.0pi*gamma)
    y = r * sin(2.0pi*gamma)
    return PVector2D(x, y)
end

"""
$(TYPEDSIGNATURES)
Generate `PVector` from radius
"""
function rand_pos_3d(r::Number)
    eta, gamma = rand(2)

    x = r * sin(acos(2.0eta-1.0)) * cos(2.0pi*gamma)
    y = r * sin(acos(2.0eta-1.0)) * sin(2.0pi*gamma)
    z = r * (2.0eta-1.0)
    return PVector(x, y, z)
end

"""
$(TYPEDSIGNATURES)
Rotate along the positive z-axis.
- `ratio`: rotational motion v.s. random motion. If equals `1`, only rotational component
"""
function rotational_velocity_acc(x, y, z, a, rot_ratio = 1.0)
    u = unit(x)
    r = sqrt(x^2 + y^2 + z^2)
    v = sqrt(a * r)
    v_vec_rot = -normalize(PVector(ustrip(u,x), ustrip(u,y), 0.0) × PVector(0.0, 0.0, 1.0)) * v
    v_vec = rot_ratio * v_vec_rot + (1-rot_ratio) * normalize(randn(PVector{Float64})) * v 
    # v_vec = v_vec_rot + (1-rot_ratio) * randn(PVector{Float64}) * v # this will change the average velocity
end

"""
$(TYPEDSIGNATURES)
Rotate along the positive z-axis.
- `rot_ratio`: rotational motion v.s. random motion. If equals `1`, only rotational component
"""
function rotational_velocity(x, y, v, rot_ratio = 1.0)
    u = unit(x)
    v_vec_rot = -normalize(PVector(ustrip(u,x), ustrip(u,y), 0.0) × PVector(0.0, 0.0, 1.0)) * v
    v_vec = rot_ratio * v_vec_rot + (1-rot_ratio) * normalize(randn(PVector{Float64})) * v 
    # v_vec = v_vec_rot + (1-rot_ratio) * randn(PVector{Float64}) * v # this will change the average velocity
end


"""
$(TYPEDSIGNATURES)
"""
function freefall_velocity_acc(x, y, z, a)
    u = unit(x)
    r = sqrt(x^2 + y^2 + z^2)
    v = sqrt(a * r)
    v_vec = -normalize(PVector(ustrip(u,x), ustrip(u,y), ustrip(u,z))) * v
end