```@meta
CurrentModule = QuantumESPRESSO
```

# QuantumESPRESSO

Documentation for [QuantumESPRESSO](https://github.com/MineralsCloud/QuantumESPRESSO.jl).

QuantumESPRESSO.jl is simply a wrapper of the types, methods, and commands defined in
[QuantumESPRESSOBase.jl](https://github.com/MineralsCloud/QuantumESPRESSOBase.jl),
[QuantumESPRESSOParser.jl](https://github.com/MineralsCloud/QuantumESPRESSOParser.jl),
[QuantumESPRESSOFormatter.jl](https://github.com/MineralsCloud/QuantumESPRESSOFormatter.jl),
and [QuantumESPRESSOCommands.jl](https://github.com/MineralsCloud/QuantumESPRESSOCommands.jl)
under a common namespace.

See the [Index](@ref main-index) for the complete list of documented functions
and types.

The code is [hosted on GitHub](https://github.com/MineralsCloud/QuantumESPRESSO.jl),
with some continuous integration services to test its validity.

This repository is created and maintained by [@singularitti](https://github.com/singularitti).
You are very welcome to contribute.

## Installation

The package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
pkg> add QuantumESPRESSO
```

Or, equivalently, via the `Pkg` API:

```@repl
import Pkg; Pkg.add("QuantumESPRESSO")
```

## Documentation

- [**STABLE**](https://MineralsCloud.github.io/QuantumESPRESSO.jl/stable) — **documentation of the most recently tagged version.**
- [**DEV**](https://MineralsCloud.github.io/QuantumESPRESSO.jl/dev) — _documentation of the in-development version._

## Project status

The package is tested against, and being developed for, Julia `1.6` and above on Linux,
macOS, and Windows.

## Questions and contributions

Usage questions can be posted on
[our discussion page](https://github.com/MineralsCloud/QuantumESPRESSO.jl/discussions).

Contributions are very welcome, as are feature requests and suggestions. Please open an
[issue](https://github.com/MineralsCloud/QuantumESPRESSO.jl/issues)
if you encounter any problems. The [Contributing](@ref) page has
a few guidelines that should be followed when opening pull requests and contributing code.

## Manual outline

```@contents
Pages = [
    "installation.md",
    "developers/contributing.md",
    "developers/style-guide.md",
    "developers/design-principles.md",
    "troubleshooting.md",
]
Depth = 3
```

## Library outline

```@contents
Pages = ["public.md"]
```

### [Index](@id main-index)

```@index
Pages = ["public.md"]
```
