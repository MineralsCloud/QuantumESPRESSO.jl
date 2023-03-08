module Outputs

using Reexport: @reexport

@reexport using QuantumESPRESSOParser.Outputs: SubroutineError

module PWscf

using Reexport: @reexport

@reexport using QuantumESPRESSOParser.Outputs.PWscf

end

module PHonon

using Reexport: @reexport

@reexport using QuantumESPRESSOParser.Outputs.PHonon

end

end
