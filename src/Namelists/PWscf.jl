"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using Parameters: @with_kw

using QuantumESPRESSO.Namelists

@with_kw struct ControlNamelist <: Namelist
    calculation::String = "scf"
    title::String = " "
    verbosity::String = "low"
    restart_mode::String = "from_scratch"
    wf_collect::Bool = true
    nstep::Int = 1
    iprint::Int = 1
    tstress::Bool = false
    tprnfor::Bool = false
    dt::Float64 = 20.0
    outdir::String = "./"
    wfcdir::String = "./"
    prefix::String = "pwscf"
    lkpoint_dir::Bool = true
    max_seconds::Float64 = 10000000.0
    etot_conv_thr::Float64 = 0.0001
    forc_conv_thr::Float64 = 0.001
    disk_io::String = "medium"
    pseudo_dir::String = "\$HOME/espresso/pseudo/"
    tefield::Bool = false
    dipfield::Bool = false
    lelfield::Bool = false
    nberrycyc::Int = 1
    lorbm::Bool = false
    lberry::Bool = false
    gdir::Int = 1
    nppstr::Int = 1
    lfcpopt::Bool = false
    gate::Bool = false
end

end
