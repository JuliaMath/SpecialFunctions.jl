using Compat
using Compat.Sys: isapple, islinux, iswindows

did_setup = false

if VERSION < v"0.7.0-DEV.1760"
    # No need to build or download anything; openspecfun is part of Julia
elseif get(ENV, "JULIA_SPECIALFUNCTIONS_BUILD_SOURCE", "false") == "true"
    # Allow a fast-path for building from source
    info("Building openspecfun from source by request")
    include("scratch.jl")
elseif isapple() || iswindows()
    # Windows and macOS can always use our binaries, and we have no binaries
    # for other non-Linux systems (e.g. BSDs)
    include("binaries.jl")
elseif !islinux()
    include("scratch.jl")
else # linux
    # Determine the glibc version. If the check fails, we know we're on a non-glibc
    # system, which means we can't use the binaries and need to build from source.
    # The glibc version used by the binaries is 2.6, so we need at least that.
    libc_ptr = ccall(:jl_dlopen, Ptr{Void}, (Ptr{Void}, UInt32), C_NULL, 0)
    glibc_ptr = Libdl.dlsym_e(libc_ptr, :gnu_get_libc_version)
    if glibc_ptr == C_NULL
        include("scratch.jl")
    else
        glibc_vers = unsafe_string(ccall(glibc_ptr, Ptr{UInt8}, ()))
        if isempty(glibc_vers) || VersionNumber(glibc_vers) < v"2.6.0"
            include("scratch.jl")
        else
            include("binaries.jl")
        end
    end
end
