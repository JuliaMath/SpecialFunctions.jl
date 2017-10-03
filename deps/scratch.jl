# Building OpenSpecFun from scratch

using BinDeps
using BinDeps: libdir, srcdir, includedir, depsdir, builddir
using Base.Math: libm

# Don't call setup again if we're being included from binaries.jl
if !did_setup
    BinDeps.@setup
    const OSF_VERS = v"0.5.3"
    openspecfun = library_dependency("libopenspecfun")
end

# If Julia is built with OpenLibm, we want to build OpenSpecFun with it as well.
# Unfortunately this requires a fair bit more work, as we need to link to the .so
# and to include the headers, which aren't readily available.
if libm == "libopenlibm"
    const OLM_VERS = v"0.5.4"
    use_openlibm = true

    if !isdir(libdir(openspecfun))
        mkpath(libdir(openspecfun))
    end

    openlibm_so = Libdl.dlpath(libm)

    # Copy over the OpenLibm .so
    cp(openlibm_so, joinpath(libdir(openspecfun), basename(openlibm_so)),
       remove_destination=true, follow_symlinks=true)

    if !isdir(srcdir(openspecfun))
        mkpath(srcdir(openspecfun))
    end

    # Grab and unpack the tarball so we can get the header files
    openlibm_tarball = joinpath(srcdir(openspecfun), "openlibm-$OLM_VERS.tar.gz")
    run(```
        curl -fkL --connect-timeout 15 -y 15
        https://github.com/JuliaLang/openlibm/archive/v$OLM_VERS.tar.gz
        -o $openlibm_tarball
    ```)
    openlibm_src = joinpath(srcdir(openspecfun), "openlibm")
    if !isdir(openlibm_src)
        mkpath(openlibm_src)
    end
    run(`tar -C $openlibm_src --strip-components 1 -xf $openlibm_tarball`)

    # Copy over all of the OpenLibm headers
    openlibm_include = joinpath(includedir(openspecfun), "openlibm")
    if !isdir(openlibm_include)
        mkpath(openlibm_include)
    end
    for f in readdir(joinpath(openlibm_src, "include"))
        cp(joinpath(openlibm_src, "include", f), joinpath(openlibm_include, f),
           remove_destination=true, follow_symlinks=true)
    end
    for f in readdir(joinpath(openlibm_src, "src"))
        if endswith(f, ".h")
            cp(joinpath(openlibm_src, "src", f), joinpath(openlibm_include, f),
               remove_destination=true, follow_symlinks=true)
        end
    end
else
    use_openlibm = false
end

fc = "gfortran"

# macOS has precompiled binaries, so it's just FreeBSD that should default to Clang
if Sys.KERNEL === :FreeBSD
    cc = "clang"
    use_clang = true
else
    cc = "gcc"
    use_clang = false
end

if Sys.ARCH in [:i386, :i387, :i486, :i586, :i686]
    cc *= " -m32"
    fc *= " -m32"
elseif Sys.ARCH === :x86_64
    cc *= " -m64"
    fc *= " -m64"
end

flags = [
    # OpenSpecFun build flags
    "ARCH=\"$(Sys.ARCH)\"",
    "CC=\"$cc\"",
    "FC=\"$fc\"",
    "USECLANG=$(Int(use_clang))",
    "USEGCC=$(Int(!use_clang))",
    "USE_OPENLIBM=$(Int(use_openlibm))",
    "CFLAGS=\"-O3 -std=c99\"",
    "FFLAGS=\"-O2 -fPIC\"",
    "LDFLAGS=\"-L$(libdir(openspecfun)) -Wl,-rpath,'\$\$ORIGIN' -Wl,-z,origin\"",
    # Make flags
    "DESTDIR=\"\"",
    "prefix=$(depsdir(openspecfun))",
    "libdir=$(libdir(openspecfun))",
    "shlibdir=$(libdir(openspecfun))",
    "includedir=$(includedir(openspecfun))",
    "O="
]

provides(Sources, URI("https://github.com/JuliaLang/openspecfun/archive/v$OSF_VERS.tar.gz"),
         openspecfun)

provides(BuildProcess,
    (@build_steps begin
        GetSources(openspecfun)
        CreateDirectory(builddir(openspecfun))
        @build_steps begin
            ChangeDirectory(builddir(openspecfun))
            FileRule(joinpath(libdir(openspecfun), "libopenspecfun." * Libdl.dlext),
                @build_steps begin
                    CreateDirectory(libdir(openspecfun))
                    `$MAKE_CMD install $flags`
                end)
        end
    end), openspecfun)

# If we're being included, the installation step happens once we return back
# to binaries.jl
if !did_setup
    BinDeps.@install Dict(:libopenspecfun => :openspecfun)
end
