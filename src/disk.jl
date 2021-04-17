mutable struct ExponentialDisk{I, Len, MASS} <: InitialConditionConfig
    collection::Collection
    NumSamples::I

    TotalMass::MASS

    ScaleRadius::Len
    ScaleHeight::Len
end

"""
    struct ExponentialDisk

## Fields

- `collection` particle type
- `NumSamples` amount of particles
"""
function ExponentialDisk(;
        collection::Collection = STAR,
        NumSamples::Int64 = 1000,
        TotalMass::Number = 1.0e10u"Msun",
        
    )
    
end

function Base.show(io::IO, config::ExponentialDisk)
    print(io,
        """
        Config of Exponential Disk:

        """
    )
end

function generate(config)
    
end