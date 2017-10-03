using BinDeps, Compat
using BinDeps: libdir

modified_defaults = false
if !in(BinDeps.Binaries, BinDeps.defaults)
    unshift!(BinDeps.defaults, BinDeps.Binaries)
    modified_defaults = true
end

BinDeps.@setup
did_setup = true

const OSF_VERS = v"0.5.3"

openspecfun = library_dependency("libopenspecfun")

const URL = "https://github.com/ararslan/openspecfun-builder/releases/download/v$OSF_VERS" *
            "/libopenspecfun-$OSF_VERS"

const DOWNLOADS = Dict(
    "Linux-x86_64"   => ("$URL-linux-x86_64.tar.gz",
                         "d70a2a391915f64f44da21915bf93ce08d054127028088addca36e16ac53bcb1"),
    "Linux-i686"     => ("$URL-linux-i686.tar.gz",
                         "e5418b170b537af2f7f1f1d06eee9be01555404f5d22a47e18bc06a540321478"),
    "Darwin-x86_64"  => ("$URL-osx-x86_64.tar.gz",
                         "e57f5f84439757a2fd1d3821a6e19a3fa69b5b1e181cc40fec0d1652fbb9efdc"),
    "Windows-x86_64" => ("$URL-win-x86_64.zip",
                         "7a5f7be4ed46d7f9d6d18a599157075512c50a372da2b2908079a3dcab9a0f25"),
    "Windows-i686"   => ("$URL-win-i686.zip",
                         "2f63a08d80e67964e2c368367f4caef7039080828e217d288669416cd46f4584"),
)

const SYSTEM = string(BinDeps.OSNAME, '-', Sys.ARCH)

if haskey(DOWNLOADS, SYSTEM)
    url, sha = DOWNLOADS[SYSTEM]
    provides(Binaries, URI(url), openspecfun, SHA=sha, os=BinDeps.OSNAME,
             unpacked_dir=joinpath("usr", "lib"), installed_libpath=libdir(openspecfun))
else
    info("No precompiled binaries found for your system. Building from scratch...")
    include("scratch.jl")
end

BinDeps.@install Dict(:libopenspecfun => :openspecfun)

if modified_defaults
    shift!(BinDeps.defaults)
end
