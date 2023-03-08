module Inputs

using Reexport: @reexport

@reexport using QuantumESPRESSOBase.Inputs
@reexport using QuantumESPRESSOParser.Inputs: InvalidInput

module PWscf

using Reexport: @reexport

@reexport using QuantumESPRESSOBase.Inputs.PWscf
@reexport using QuantumESPRESSOParser.Inputs.PWscf
@reexport using QuantumESPRESSOFormatter.Inputs.PWscf

end

module PHonon

using Reexport: @reexport

@reexport using QuantumESPRESSOBase.Inputs.PHonon
@reexport using QuantumESPRESSOParser.Inputs.PHonon

end

end
