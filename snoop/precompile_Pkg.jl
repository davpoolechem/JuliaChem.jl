function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    isdefined(Pkg, Symbol("##ensure_artifact_installed#42")) && precompile(Tuple{getfield(Pkg.Artifacts, Symbol("##ensure_artifact_installed#42")), Pkg.BinaryPlatforms.Platform, Bool, Bool, typeof(Pkg.Artifacts.ensure_artifact_installed), String, Base.Dict{String, Any}, String})
    isdefined(Pkg, Symbol("##ensure_artifact_installed#42")) && precompile(Tuple{getfield(Pkg.Artifacts, Symbol("##ensure_artifact_installed#42")), Pkg.BinaryPlatforms.Platform, Bool, Bool, typeof(Pkg.Artifacts.ensure_artifact_installed), String, Base.Dict{String, Any}, String})
    isdefined(Pkg, Symbol("##query_override#7")) && precompile(Tuple{getfield(Pkg.Artifacts, Symbol("##query_override#7")), Base.Dict{Symbol, Base.Dict{K, V} where V where K}, typeof(Pkg.Artifacts.query_override), Base.SHA1})
    isdefined(Pkg, Symbol("##query_override#7")) && precompile(Tuple{getfield(Pkg.Artifacts, Symbol("##query_override#7")), Base.Dict{Symbol, Base.Dict{K, V} where V where K}, typeof(Pkg.Artifacts.query_override), Base.SHA1})
    isdefined(Pkg, Symbol("#ensure_artifact_installed##kw")) && precompile(Tuple{getfield(Pkg.Artifacts, Symbol("#ensure_artifact_installed##kw")), NamedTuple{(:platform,), Tuple{Pkg.BinaryPlatforms.Linux}}, typeof(Pkg.Artifacts.ensure_artifact_installed), String, Base.Dict{String, Any}, String})
    isdefined(Pkg, Symbol("#ensure_artifact_installed##kw")) && precompile(Tuple{getfield(Pkg.Artifacts, Symbol("#ensure_artifact_installed##kw")), NamedTuple{(:platform,), Tuple{Pkg.BinaryPlatforms.Linux}}, typeof(Pkg.Artifacts.ensure_artifact_installed), String, Base.Dict{String, Any}, String})
    precompile(Tuple{typeof(Pkg.Artifacts.do_artifact_str), String, Base.Dict{String, Any}, String, Module})
    precompile(Tuple{typeof(Pkg.Artifacts.do_artifact_str), String, Base.Dict{String, Any}, String, Module})
end
