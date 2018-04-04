var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Special-Functions-1",
    "page": "Home",
    "title": "Special Functions",
    "category": "section",
    "text": "This package provides a comprehensive collection of special functions based on the OpenSpecFun and OpenLibm libraries.Function Description\nerf(x) error function at x\nerfc(x) complementary error function, i.e. the accurate version of 1-erf(x) for large x\nerfinv(x) inverse function to erf()\nerfcinv(x) inverse function to erfc()\nerfi(x) imaginary error function defined as -im * erf(x * im), where im is the imaginary unit\nerfcx(x) scaled complementary error function, i.e. accurate exp(x^2) * erfc(x) for large x\ndawson(x) scaled imaginary error function, a.k.a. Dawson function, i.e. accurate exp(-x^2) * erfi(x) * sqrt(pi) / 2 for large x\nsinint(x) sine integral at x\ncosint(x) cosine integral at x\ndigamma(x) digamma function (i.e. the derivative of lgamma at x)\neta(x) Dirichlet eta function at x\nzeta(x) Riemann zeta function at x\nairyai(z) Airy Ai function at z\nairyaiprime(z) derivative of the Airy Ai function at z\nairybi(z) Airy Bi function at z\nairybiprime(z) derivative of the Airy Bi function at z\nairyaix(z), airyaiprimex(z), airybix(z), airybiprimex(z) scaled Airy Ai function and kth derivatives at z\nbesselj(nu,z) Bessel function of the first kind of order nu at z\nbesselj0(z) besselj(0,z)\nbesselj1(z) besselj(1,z)\nbesseljx(nu,z) scaled Bessel function of the first kind of order nu at z\nbessely(nu,z) Bessel function of the second kind of order nu at z\nbessely0(z) bessely(0,z)\nbessely1(z) bessely(1,z)\nbesselyx(nu,z) scaled Bessel function of the second kind of order nu at z\nbesselh(nu,k,z) Bessel function of the third kind (a.k.a. Hankel function) of order nu at z; k must be either 1 or 2\nhankelh1(nu,z) besselh(nu, 1, z)\nhankelh1x(nu,z) scaled besselh(nu, 1, z)\nhankelh2(nu,z) besselh(nu, 2, z)\nhankelh2x(nu,z) scaled besselh(nu, 2, z)\nbesseli(nu,z) modified Bessel function of the first kind of order nu at z\nbesselix(nu,z) scaled modified Bessel function of the first kind of order nu at z\nbesselk(nu,z) modified Bessel function of the second kind of order nu at z\nbesselkx(nu,z) scaled modified Bessel function of the second kind of order nu at z"
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "The package is available for Julia versions 0.5 and up. To install it, runPkg.add(\"SpecialFunctions\")from the Julia REPL."
},

{
    "location": "index.html#Note-1",
    "page": "Home",
    "title": "Note",
    "category": "section",
    "text": "Prior to Julia 0.6, most of these functions were available in Julia\'s Base module. Because of this, the symbols from this package are not exported on Julia 0.5 to avoid name conflicts. In this case, the symbols will need to be explicitly imported or called with the prefix SpecialFunctions. This is not necessary for Julia versions 0.6 and later.On Julia 0.7, openspecfun is not build as part of Julia. Thus for Julia versions 0.7 and later, installing this package downloads openspecfun. Binaries of openspecfun are available for macOS, Windows, and Linux (glibc >= 2.6). Other systems will need to build the library from source. You can force a build from source by setting an environment variable called JULIA_SPECIALFUNCTIONS_BUILD_SOURCE equal to true before running Pkg.build. This ensures that the library is built locally from source, even if binaries are available. Doing this requires a C compiler (Clang on macOS and FreeBSD, GCC elsewhere) and gfortran. If you always want to build this library from source, consider addingENV[\"JULIA_SPECIALFUNCTIONS_BUILD_SOURCE\"] = \"true\"to your .juliarc.jl file."
},

{
    "location": "special.html#",
    "page": "Functions",
    "title": "Functions",
    "category": "page",
    "text": ""
},

{
    "location": "special.html#SpecialFunctions.erf",
    "page": "Functions",
    "title": "SpecialFunctions.erf",
    "category": "function",
    "text": "erf(x)\n\nCompute the error function of x, defined by frac2sqrtpi int_0^x e^-t^2 dt for arbitrary complex x.\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.erfc",
    "page": "Functions",
    "title": "SpecialFunctions.erfc",
    "category": "function",
    "text": "erfc(x)\n\nCompute the complementary error function of x, defined by 1 - operatornameerf(x).\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.erfcx",
    "page": "Functions",
    "title": "SpecialFunctions.erfcx",
    "category": "function",
    "text": "erfcx(x)\n\nCompute the scaled complementary error function of x, defined by e^x^2 operatornameerfc(x). Note also that operatornameerfcx(-ix) computes the Faddeeva function w(x).\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.erfi",
    "page": "Functions",
    "title": "SpecialFunctions.erfi",
    "category": "function",
    "text": "erfi(x)\n\nCompute the imaginary error function of x, defined by -i operatornameerf(ix).\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.dawson",
    "page": "Functions",
    "title": "SpecialFunctions.dawson",
    "category": "function",
    "text": "dawson(x)\n\nCompute the Dawson function (scaled imaginary error function) of x, defined by fracsqrtpi2 e^-x^2 operatornameerfi(x).\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.erfinv",
    "page": "Functions",
    "title": "SpecialFunctions.erfinv",
    "category": "function",
    "text": "erfinv(x)\n\nCompute the inverse error function of a real x, defined by operatornameerf(operatornameerfinv(x)) = x.\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.erfcinv",
    "page": "Functions",
    "title": "SpecialFunctions.erfcinv",
    "category": "function",
    "text": "erfcinv(x)\n\nCompute the inverse error complementary function of a real x, defined by operatornameerfc(operatornameerfcinv(x)) = x.\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.sinint",
    "page": "Functions",
    "title": "SpecialFunctions.sinint",
    "category": "function",
    "text": "sinint(x)\n\nCompute the sine integral function of x, defined by operatornameSi(x) = int_0^xfracsin tt dt for real x.\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.cosint",
    "page": "Functions",
    "title": "SpecialFunctions.cosint",
    "category": "function",
    "text": "cosint(x)\n\nCompute the cosine integral function of x, defined by operatornameCi(x) = gamma + log x + int_0^x fraccos t - 1t dt for real x > 0, where gamma is the Euler-Mascheroni constant.\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.digamma",
    "page": "Functions",
    "title": "SpecialFunctions.digamma",
    "category": "function",
    "text": "digamma(x)\n\nCompute the digamma function of x (the logarithmic derivative of gamma(x)).\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.invdigamma",
    "page": "Functions",
    "title": "SpecialFunctions.invdigamma",
    "category": "function",
    "text": "invdigamma(x)\n\nCompute the inverse digamma function of x.\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.trigamma",
    "page": "Functions",
    "title": "SpecialFunctions.trigamma",
    "category": "function",
    "text": "trigamma(x)\n\nCompute the trigamma function of x (the logarithmic second derivative of gamma(x)).\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.polygamma",
    "page": "Functions",
    "title": "SpecialFunctions.polygamma",
    "category": "function",
    "text": "polygamma(m, x)\n\nCompute the polygamma function of order m of argument x (the (m+1)th derivative of the logarithm of gamma(x))\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.airyai",
    "page": "Functions",
    "title": "SpecialFunctions.airyai",
    "category": "function",
    "text": "airyai(x)\n\nAiry function of the first kind operatornameAi(x).\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.airyaiprime",
    "page": "Functions",
    "title": "SpecialFunctions.airyaiprime",
    "category": "function",
    "text": "airyaiprime(x)\n\nDerivative of the Airy function of the first kind operatornameAi(x).\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.airyaix",
    "page": "Functions",
    "title": "SpecialFunctions.airyaix",
    "category": "function",
    "text": "airyaix(x)\n\nScaled Airy function of the first kind operatornameAi(x) e^frac23 x sqrtx.  Throws DomainError for negative Real arguments.\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.airyaiprimex",
    "page": "Functions",
    "title": "SpecialFunctions.airyaiprimex",
    "category": "function",
    "text": "airyaiprimex(x)\n\nScaled derivative of the Airy function of the first kind operatornameAi(x) e^frac23 x sqrtx.  Throws DomainError for negative Real arguments.\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.airybi",
    "page": "Functions",
    "title": "SpecialFunctions.airybi",
    "category": "function",
    "text": "airybi(x)\n\nAiry function of the second kind operatornameBi(x).\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.airybiprime",
    "page": "Functions",
    "title": "SpecialFunctions.airybiprime",
    "category": "function",
    "text": "airybiprime(x)\n\nDerivative of the Airy function of the second kind operatornameBi(x).\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.airybix",
    "page": "Functions",
    "title": "SpecialFunctions.airybix",
    "category": "function",
    "text": "airybix(x)\n\nScaled Airy function of the second kind operatornameBi(x) e^- left operatornameRe left( frac23 x sqrtx right) right.\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.airybiprimex",
    "page": "Functions",
    "title": "SpecialFunctions.airybiprimex",
    "category": "function",
    "text": "airybiprimex(x)\n\nScaled derivative of the Airy function of the second kind operatornameBi(x) e^- left operatornameRe left( frac23 x sqrtx right) right.\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.besselj0",
    "page": "Functions",
    "title": "SpecialFunctions.besselj0",
    "category": "function",
    "text": "besselj0(x)\n\nBessel function of the first kind of order 0, J_0(x).\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.besselj1",
    "page": "Functions",
    "title": "SpecialFunctions.besselj1",
    "category": "function",
    "text": "besselj1(x)\n\nBessel function of the first kind of order 1, J_1(x).\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.besselj",
    "page": "Functions",
    "title": "SpecialFunctions.besselj",
    "category": "function",
    "text": "besselj(nu, x)\n\nBessel function of the first kind of order nu, J_nu(x).\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.besseljx",
    "page": "Functions",
    "title": "SpecialFunctions.besseljx",
    "category": "function",
    "text": "besseljx(nu, x)\n\nScaled Bessel function of the first kind of order nu, J_nu(x) e^-  operatornameIm(x) .\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.bessely0",
    "page": "Functions",
    "title": "SpecialFunctions.bessely0",
    "category": "function",
    "text": "bessely0(x)\n\nBessel function of the second kind of order 0, Y_0(x).\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.bessely1",
    "page": "Functions",
    "title": "SpecialFunctions.bessely1",
    "category": "function",
    "text": "bessely1(x)\n\nBessel function of the second kind of order 1, Y_1(x).\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.bessely",
    "page": "Functions",
    "title": "SpecialFunctions.bessely",
    "category": "function",
    "text": "bessely(nu, x)\n\nBessel function of the second kind of order nu, Y_nu(x).\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.besselyx",
    "page": "Functions",
    "title": "SpecialFunctions.besselyx",
    "category": "function",
    "text": "besselyx(nu, x)\n\nScaled Bessel function of the second kind of order nu, Y_nu(x) e^-  operatornameIm(x) .\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.hankelh1",
    "page": "Functions",
    "title": "SpecialFunctions.hankelh1",
    "category": "function",
    "text": "hankelh1(nu, x)\n\nBessel function of the third kind of order nu, H^(1)_nu(x).\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.hankelh1x",
    "page": "Functions",
    "title": "SpecialFunctions.hankelh1x",
    "category": "function",
    "text": "hankelh1x(nu, x)\n\nScaled Bessel function of the third kind of order nu, H^(1)_nu(x) e^-x i.\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.hankelh2",
    "page": "Functions",
    "title": "SpecialFunctions.hankelh2",
    "category": "function",
    "text": "hankelh2(nu, x)\n\nBessel function of the third kind of order nu, H^(2)_nu(x).\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.hankelh2x",
    "page": "Functions",
    "title": "SpecialFunctions.hankelh2x",
    "category": "function",
    "text": "hankelh2x(nu, x)\n\nScaled Bessel function of the third kind of order nu, H^(2)_nu(x) e^x i.\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.besselh",
    "page": "Functions",
    "title": "SpecialFunctions.besselh",
    "category": "function",
    "text": "besselh(nu, [k=1,] x)\n\nBessel function of the third kind of order nu (the Hankel function). k is either 1 or 2, selecting hankelh1 or hankelh2, respectively. k defaults to 1 if it is omitted. (See also besselhx for an exponentially scaled variant.)\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.besselhx",
    "page": "Functions",
    "title": "SpecialFunctions.besselhx",
    "category": "function",
    "text": "besselhx(nu, [k=1,] z)\n\nCompute the scaled Hankel function exp(iz) H_^(k)(z), where k is 1 or 2, H_^(k)(z) is besselh(nu, k, z), and  is - for k=1 and + for k=2.  k defaults to 1 if it is omitted.\n\nThe reason for this function is that H_^(k)(z) is asymptotically proportional to exp(iz)sqrtz for large z, and so the besselh function is susceptible to overflow or underflow when z has a large imaginary part.  The besselhx function cancels this exponential factor (analytically), so it avoids these problems.\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.besseli",
    "page": "Functions",
    "title": "SpecialFunctions.besseli",
    "category": "function",
    "text": "besseli(nu, x)\n\nModified Bessel function of the first kind of order nu, I_nu(x).\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.besselix",
    "page": "Functions",
    "title": "SpecialFunctions.besselix",
    "category": "function",
    "text": "besselix(nu, x)\n\nScaled modified Bessel function of the first kind of order nu, I_nu(x) e^-  operatornameRe(x) .\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.besselk",
    "page": "Functions",
    "title": "SpecialFunctions.besselk",
    "category": "function",
    "text": "besselk(nu, x)\n\nModified Bessel function of the second kind of order nu, K_nu(x).\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.besselkx",
    "page": "Functions",
    "title": "SpecialFunctions.besselkx",
    "category": "function",
    "text": "besselkx(nu, x)\n\nScaled modified Bessel function of the second kind of order nu, K_nu(x) e^x.\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.eta",
    "page": "Functions",
    "title": "SpecialFunctions.eta",
    "category": "function",
    "text": "eta(x)\n\nDirichlet eta function eta(s) = sum^infty_n=1(-1)^n-1n^s.\n\n\n\n\n\n"
},

{
    "location": "special.html#SpecialFunctions.zeta",
    "page": "Functions",
    "title": "SpecialFunctions.zeta",
    "category": "function",
    "text": "zeta(s, z)\n\nGeneralized zeta function zeta(s z), defined by the sum sum_k=0^infty ((k+z)^2)^-s2, where any term with k+z=0 is excluded.  For Re z  0, this definition is equivalent to the Hurwitz zeta function sum_k=0^infty (k+z)^-s.   For z=1, it yields the Riemann zeta function zeta(s).\n\n\n\n\n\nzeta(s)\n\nRiemann zeta function zeta(s).\n\n\n\n\n\n"
},

{
    "location": "special.html#Functions-1",
    "page": "Functions",
    "title": "Functions",
    "category": "section",
    "text": "CurrentModule = SpecialFunctionsSpecialFunctions.erf\nSpecialFunctions.erfc\nSpecialFunctions.erfcx\nSpecialFunctions.erfi\nSpecialFunctions.dawson\nSpecialFunctions.erfinv\nSpecialFunctions.erfcinv\nSpecialFunctions.sinint\nSpecialFunctions.cosint\nSpecialFunctions.digamma\nSpecialFunctions.invdigamma\nSpecialFunctions.trigamma\nSpecialFunctions.polygamma\nSpecialFunctions.airyai\nSpecialFunctions.airyaiprime\nSpecialFunctions.airyaix\nSpecialFunctions.airyaiprimex\nSpecialFunctions.airybi\nSpecialFunctions.airybiprime\nSpecialFunctions.airybix\nSpecialFunctions.airybiprimex\nSpecialFunctions.besselj0\nSpecialFunctions.besselj1\nSpecialFunctions.besselj\nSpecialFunctions.besseljx\nSpecialFunctions.bessely0\nSpecialFunctions.bessely1\nSpecialFunctions.bessely\nSpecialFunctions.besselyx\nSpecialFunctions.hankelh1\nSpecialFunctions.hankelh1x\nSpecialFunctions.hankelh2\nSpecialFunctions.hankelh2x\nSpecialFunctions.besselh\nSpecialFunctions.besselhx\nSpecialFunctions.besseli\nSpecialFunctions.besselix\nSpecialFunctions.besselk\nSpecialFunctions.besselkx\nSpecialFunctions.eta\nSpecialFunctions.zeta"
},

]}
