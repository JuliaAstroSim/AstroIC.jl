# Tools

Here are some tools to manipulate initial conditions

## Example: collisions

```julia
config = PlummerStarCluster(NumSamples = 1000)
galaxy1 = generate(config);
galaxy2 = generate(config);
galaxy3 = generate(config);
galaxy4 = generate(config);

setpos(galaxy1, PVector(-0.1, 0.0, 0.0, u"kpc"))
setpos(galaxy2, PVector(+0.1, 0.0, 0.0, u"kpc"))
setpos(galaxy3, PVector(0.0, +0.1, 0.0, u"kpc"))
setpos(galaxy4, PVector(0.0, -0.1, 0.0, u"kpc"))

setvel(galaxy1, PVector(+0.4, +0.8, 0.0, u"kpc/Gyr"))
setvel(galaxy2, PVector(-0.4, -1.0, -0.2, u"kpc/Gyr"))
setvel(galaxy3, PVector(+1.0, -0.4, 0.4, u"kpc/Gyr"))
setvel(galaxy4, PVector(-1.0, +0.4, -0.2, u"kpc/Gyr"))

data = deepcopy(galaxy1);
append!(data, galaxy2);
append!(data, galaxy3);
append!(data, galaxy4);
data.ID .= 1:4000;
```