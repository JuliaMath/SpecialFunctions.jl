using BinDeps

BinDeps.@setup

const OSF_VERSION = v"0.5.3"
const SOMAJOR, SOMINOR = 1, 3 # NOTE: Ensure these match JuliaLang/openspecfun/Make.inc

openspecfun = library_dependency("openspecfun", aliases=["libopenspecfun.$SOMAJOR", "libopenspecfun.$SOMAJOR.$SOMINOR",
                                                         "libopenspecfun.so.$SOMAJOR", "libopenspecfun.so.$SOMAJOR.$SOMINOR",
                                                         "libopenspecfun"])

prefix = BinDeps.usrdir(openspecfun)
source = joinpath(BinDeps.srcdir(openspecfun), "openspecfun-$OSF_VERSION")
shlib = joinpath(BinDeps.libdir(openspecfun), "libopenspecfun." * BinDeps.shlib_ext)
make = is_bsd() && !is_apple() ? "gmake" : "make"

provides(Sources, URI("https://github.com/JuliaLang/openspecfun/archive/v$OSF_VERSION.tar.gz"),
         openspecfun, unpacked_dir="openspecfun-$OSF_VERSION")

if is_apple()
    using Homebrew
    provides(Homebrew.HB, "staticfloat/juliadeps/openspecfun", openspecfun, os=:Darwin)
end

provides(BuildProcess, (@build_steps begin
    GetSources(openspecfun)
    CreateDirectory(BinDeps.libdir(openspecfun))
    @build_steps begin
        ChangeDirectory(source)
        FileRule(shlib, @build_steps begin
            `$make install prefix=$prefix`
        end)
    end
end), openspecfun)

BinDeps.@install Dict(:openspecfun => :openspecfun)
