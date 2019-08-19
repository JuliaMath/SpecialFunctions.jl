var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#Special-Functions-1",
    "page": "Home",
    "title": "Special Functions",
    "category": "section",
    "text": "This package provides a comprehensive collection of special functions based on the OpenSpecFun and OpenLibm libraries."
},

{
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "The package is available for Julia versions 0.5 and up. To install it, runPkg.add(\"SpecialFunctions\")from the Julia REPL."
},

{
    "location": "#Note-1",
    "page": "Home",
    "title": "Note",
    "category": "section",
    "text": "Prior to Julia 0.6, most of these functions were available in Julia\'s Base module. Because of this, the symbols from this package are not exported on Julia 0.5 to avoid name conflicts. In this case, the symbols will need to be explicitly imported or called with the prefix SpecialFunctions. This is not necessary for Julia versions 0.6 and later.On Julia 0.7, openspecfun is not build as part of Julia. Thus for Julia versions 0.7 and later, installing this package downloads openspecfun. Binaries of openspecfun are available for macOS, Windows, and Linux (glibc >= 2.6). Other systems will need to build the library from source. You can force a build from source by setting an environment variable called JULIA_SPECIALFUNCTIONS_BUILD_SOURCE equal to true before running Pkg.build. This ensures that the library is built locally from source, even if binaries are available. Doing this requires a C compiler (Clang on macOS and FreeBSD, GCC elsewhere) and gfortran. If you always want to build this library from source, consider addingENV[\"JULIA_SPECIALFUNCTIONS_BUILD_SOURCE\"] = \"true\"to your .juliarc.jl file."
},

{
    "location": "functions_overview/#",
    "page": "Overview",
    "title": "Overview",
    "category": "page",
    "text": ""
},

{
    "location": "functions_overview/#Functions-1",
    "page": "Overview",
    "title": "Functions",
    "category": "section",
    "text": "Here the Special Functions are listed according to the structure of NIST Digital Library of Mathematical Functions."
},

{
    "location": "functions_overview/#[Gamma-Function](https://dlmf.nist.gov/5)-1",
    "page": "Overview",
    "title": "Gamma Function",
    "category": "section",
    "text": "Function Description\ngamma(z) gamma function Gamma(z)\ndigamma(x) digamma function (i.e. the derivative of lgamma at x)\ninvdigamma(x) invdigamma function (i.e. inverse of digamma function at x using fixed-point iteration algorithm)\ntrigamma(x) trigamma function (i.e the logarithmic second derivative of gamma at x)\npolygamma(m,x) polygamma function (i.e the (m+1)-th derivative of the lgamma function at x)\ngamma_inc(a,x,IND) incomplete gamma function ratio P(a,x) and Q(a,x) (i.e evaluates P(a,x) and Q(a,x)for accuracy specified by IND and returns tuple (p,q))\nbeta_inc(a,b,x,y) incomplete beta function ratio Ix(a,b) and Iy(a,b) (i.e evaluates Ix(a,b) and Iy(a,b) and returns tuple (p,q))\ngamma_inc_inv(a,p,q) inverse of incomplete gamma function ratio P(a,x) and Q(a,x) (i.e evaluates x given P(a,x)=p and Q(a,x)=q\nloggamma(x) accurate log(gamma(x)) for large x\nlogabsgamma(x) accurate log(abs(gamma(x))) for large x\nlgamma(x) accurate log(gamma(x)) for large x\nlogfactorial(x) accurate log(factorial(x)) for large x; same as lgamma(x+1) for x > 1, zero otherwise\nbeta(x,y) beta function at x,y\nlogbeta(x,y) accurate log(beta(x,y)) for large x or y\nlogabsbeta(x,y) accurate log(abs(beta(x,y))) for large x or y\nlogabsbinomial(x,y) accurate log(abs(beta(x,y))) for large x or y"
},

{
    "location": "functions_overview/#[Trigonometric-Integrals](https://dlmf.nist.gov/6)-1",
    "page": "Overview",
    "title": "Trigonometric Integrals",
    "category": "section",
    "text": "Function Description\nsinint(x) sine integral Si(x)\ncosint(x) cosine integral Ci(x)"
},

{
    "location": "functions_overview/#[Error-Functions,-Dawson’s-and-Fresnel-Integrals](https://dlmf.nist.gov/7)-1",
    "page": "Overview",
    "title": "Error Functions, Dawson’s and Fresnel Integrals",
    "category": "section",
    "text": "Function Description\nerf(x) error function at x\nerfc(x) complementary error function, i.e. the accurate version of 1-operatornameerf(x) for large x\nerfcinv(x) inverse function to erfc()\nerfcx(x) scaled complementary error function, i.e. accurate e^x^2 operatornameerfc(x) for large x\nerfi(x) imaginary error function defined as -i operatornameerf(ix)\nerfinv(x) inverse function to erf()\ndawson(x) scaled imaginary error function, a.k.a. Dawson function, i.e. accurate fracsqrtpi2 e^-x^2 operatornameerfi(x) for large x"
},

{
    "location": "functions_overview/#[Airy-and-Related-Functions](https://dlmf.nist.gov/9)-1",
    "page": "Overview",
    "title": "Airy and Related Functions",
    "category": "section",
    "text": "Function Description\nairyai(z) Airy Ai function at z\nairyaiprime(z) derivative of the Airy Ai function at z\nairybi(z) Airy Bi function at z\nairybiprime(z) derivative of the Airy Bi function at z\nairyaix(z), airyaiprimex(z), airybix(z), airybiprimex(z) scaled Airy Ai function and kth derivatives at z"
},

{
    "location": "functions_overview/#[Bessel-Functions](https://dlmf.nist.gov/10)-1",
    "page": "Overview",
    "title": "Bessel Functions",
    "category": "section",
    "text": "Function Description\nbesselj(nu,z) Bessel function of the first kind of order nu at z\nbesselj0(z) besselj(0,z)\nbesselj1(z) besselj(1,z)\nbesseljx(nu,z) scaled Bessel function of the first kind of order nu at z\nbessely(nu,z) Bessel function of the second kind of order nu at z\nbessely0(z) bessely(0,z)\nbessely1(z) bessely(1,z)\nbesselyx(nu,z) scaled Bessel function of the second kind of order nu at z\nbesselh(nu,k,z) Bessel function of the third kind (a.k.a. Hankel function) of order nu at z; k must be either 1 or 2\nhankelh1(nu,z) besselh(nu, 1, z)\nhankelh1x(nu,z) scaled besselh(nu, 1, z)\nhankelh2(nu,z) besselh(nu, 2, z)\nhankelh2x(nu,z) scaled besselh(nu, 2, z)\nbesseli(nu,z) modified Bessel function of the first kind of order nu at z\nbesselix(nu,z) scaled modified Bessel function of the first kind of order nu at z\nbesselk(nu,z) modified Bessel function of the second kind of order nu at z\nbesselkx(nu,z) scaled modified Bessel function of the second kind of order nu at z"
},

{
    "location": "functions_overview/#[Elliptic-Integrals](https://dlmf.nist.gov/19)-1",
    "page": "Overview",
    "title": "Elliptic Integrals",
    "category": "section",
    "text": "Function Description\nellipk(m) complete elliptic integral of 1st kind K(m)\nellipe(m) complete elliptic integral of 2nd kind E(m)"
},

{
    "location": "functions_overview/#[Zeta-and-Related-Functions](https://dlmf.nist.gov/25)-1",
    "page": "Overview",
    "title": "Zeta and Related Functions",
    "category": "section",
    "text": "Function Description\neta(x) Dirichlet eta function at x\nzeta(x) Riemann zeta function at x"
},

{
    "location": "functions_list/#",
    "page": "List",
    "title": "List",
    "category": "page",
    "text": ""
},

{
    "location": "functions_list/#SpecialFunctions.erf",
    "page": "List",
    "title": "SpecialFunctions.erf",
    "category": "function",
    "text": "erf(x)\n\nCompute the error function of x, defined by\n\noperatornameerf(x) = frac2pi int_0^x exp(-t^2)  mathrmdt\nquad textfor quad x in mathbbC  \n\nExternal links: DLMF, Wikipedia.\n\nSee also: erfc(x), erfcx(x), erfi(x), dawson(x), erfinv(x), erfcinv(x).\n\nImplementation by\n\nFloat32/Float64: C standard math library   libm.\nBigFloat: C library for multiple-precision floating-point MPFR\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.erfc",
    "page": "List",
    "title": "SpecialFunctions.erfc",
    "category": "function",
    "text": "erfc(x)\n\nCompute the complementary error function of x, defined by\n\noperatornameerfc(x)\n= 1 - operatornameerf(x)\n= frac2pi int_x^infty exp(-t^2)  mathrmdt\nquad textfor quad x in mathbbC  \n\nThis is the accurate version of 1-erf(x) for large x.\n\nExternal links: DLMF, Wikipedia.\n\nSee also: erf(x).\n\nImplementation by\n\nFloat32/Float64: C standard math library   libm.\nBigFloat: C library for multiple-precision floating-point MPFR\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.erfcx",
    "page": "List",
    "title": "SpecialFunctions.erfcx",
    "category": "function",
    "text": "erfcx(x)\n\nCompute the scaled complementary error function of x, defined by\n\noperatornameerfcx(x)\n= e^x^2 operatornameerfc(x)\nquad textfor quad x in mathbbC  \n\nThis is the accurate version of e^x^2 operatornameerfc(x) for large x. Note also that operatornameerfcx(-ix) computes the Faddeeva function w(x).\n\nExternal links: DLMF, Wikipedia.\n\nSee also: erfc(x).\n\nImplementation by\n\nFloat32/Float64: C standard math library   libm.\nBigFloat: MPFR has an open TODO item for this function until then, we use   DLMF 7.12.1 for the tail.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.erfi",
    "page": "List",
    "title": "SpecialFunctions.erfi",
    "category": "function",
    "text": "erfi(x)\n\nCompute the imaginary error function of x, defined by\n\noperatornameerfi(x)\n= -i operatornameerf(ix)\nquad textfor quad x in mathbbC  \n\nExternal links: Wikipedia.\n\nSee also: erf(x).\n\nImplementation by\n\nFloat32/Float64: C standard math library   libm.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.dawson",
    "page": "List",
    "title": "SpecialFunctions.dawson",
    "category": "function",
    "text": "dawson(x)\n\nCompute the Dawson function (scaled imaginary error function) of x, defined by\n\noperatornamedawson(x)\n= fracsqrtpi2 e^-x^2 operatornameerfi(x)\nquad textfor quad x in mathbbC  \n\nThis is the accurate version of fracsqrtpi2 e^-x^2 operatornameerfi(x) for large x.\n\nExternal links: DLMF, Wikipedia.\n\nSee also: erfi(x).\n\nImplementation by\n\nFloat32/Float64: C standard math library   libm.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.erfinv",
    "page": "List",
    "title": "SpecialFunctions.erfinv",
    "category": "function",
    "text": "erfinv(x)\n\nCompute the inverse error function of a real x, that is\n\noperatornameerfinv(x) = operatornameerf^-1(x)\nquad textfor quad x in mathbbR  \n\nExternal links: Wikipedia.\n\nSee also: erf(x).\n\nImplementation\n\nUsing the rational approximants tabulated in:\n\nJ. M. Blair, C. A. Edwards, and J. H. Johnson, \"Rational Chebyshev approximations for the inverse of the error function\", Math. Comp. 30, pp. 827–830 (1976). http://dx.doi.org/10.1090/S0025-5718-1976-0421040-7, http://www.jstor.org/stable/2005402\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.erfcinv",
    "page": "List",
    "title": "SpecialFunctions.erfcinv",
    "category": "function",
    "text": "erfcinv(x)\n\nCompute the inverse error complementary function of a real x, that is\n\noperatornameerfcinv(x) = operatornameerfc^-1(x)\nquad textfor quad x in mathbbR  \n\nExternal links: Wikipedia.\n\nSee also: erfc(x).\n\nImplementation\n\nUsing the rational approximants tabulated in:\n\nJ. M. Blair, C. A. Edwards, and J. H. Johnson, \"Rational Chebyshev approximations for the inverse of the error function\", Math. Comp. 30, pp. 827–830 (1976). http://dx.doi.org/10.1090/S0025-5718-1976-0421040-7, http://www.jstor.org/stable/2005402\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.sinint",
    "page": "List",
    "title": "SpecialFunctions.sinint",
    "category": "function",
    "text": "sinint(x)\n\nCompute the sine integral function of x, defined by\n\noperatornameSi(x)\n= int_0^x fracsin tt  mathrmdt\nquad textfor quad\nx in mathbbR \n\nExternal links: DLMF, Wikipedia.\n\nSee also: cosint(x).\n\nImplementation\n\nUsing the rational approximants tabulated in:\n\nA.J. MacLeod, \"Rational approximations, software and test methods for sine and cosine integrals\", Numer. Algor. 12, pp. 259–272 (1996). http://dx.doi.org/10.1007/BF02142806, https://link.springer.com/article/10.1007/BF02142806.\n\nNote: the second zero of Ci(x) has a typo that is fixed: r_1 = 338418 0422mathbf8 51186 42639 78511 46402 in the article, but is in fact: r_1 = 338418 0422mathbf5 51186 42639 78511 46402.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.cosint",
    "page": "List",
    "title": "SpecialFunctions.cosint",
    "category": "function",
    "text": "cosint(x)\n\nCompute the cosine integral function of x, defined by\n\noperatornameCi(x)\n= gamma + log x + int_0^x fraccos (t) - 1t  mathrmdt\nquad textfor quad\nx  0 \n\nwhere gamma is the Euler-Mascheroni constant.\n\nExternal links: DLMF, Wikipedia.\n\nSee also: sinint(x).\n\nImplementation\n\nUsing the rational approximants tabulated in:\n\nA.J. MacLeod, \"Rational approximations, software and test methods for sine and cosine integrals\", Numer. Algor. 12, pp. 259–272 (1996). http://dx.doi.org/10.1007/BF02142806, https://link.springer.com/article/10.1007/BF02142806.\n\nNote: the second zero of Ci(x) has a typo that is fixed: r_1 = 338418 0422mathbf8 51186 42639 78511 46402 in the article, but is in fact: r_1 = 338418 0422mathbf5 51186 42639 78511 46402.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.digamma",
    "page": "List",
    "title": "SpecialFunctions.digamma",
    "category": "function",
    "text": "digamma(x)\n\nCompute the digamma function of x (the logarithmic derivative of gamma(x)).\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.invdigamma",
    "page": "List",
    "title": "SpecialFunctions.invdigamma",
    "category": "function",
    "text": "invdigamma(x)\n\nCompute the inverse digamma function of x.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.trigamma",
    "page": "List",
    "title": "SpecialFunctions.trigamma",
    "category": "function",
    "text": "trigamma(x)\n\nCompute the trigamma function of x (the logarithmic second derivative of gamma(x)).\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.polygamma",
    "page": "List",
    "title": "SpecialFunctions.polygamma",
    "category": "function",
    "text": "polygamma(m, x)\n\nCompute the polygamma function of order m of argument x (the (m+1)th derivative of the logarithm of gamma(x))\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.airyai",
    "page": "List",
    "title": "SpecialFunctions.airyai",
    "category": "function",
    "text": "airyai(x)\n\nAiry function of the first kind operatornameAi(x).\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.airyaiprime",
    "page": "List",
    "title": "SpecialFunctions.airyaiprime",
    "category": "function",
    "text": "airyaiprime(x)\n\nDerivative of the Airy function of the first kind operatornameAi(x).\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.airyaix",
    "page": "List",
    "title": "SpecialFunctions.airyaix",
    "category": "function",
    "text": "airyaix(x)\n\nScaled Airy function of the first kind operatornameAi(x) e^frac23 x sqrtx.  Throws DomainError for negative Real arguments.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.airyaiprimex",
    "page": "List",
    "title": "SpecialFunctions.airyaiprimex",
    "category": "function",
    "text": "airyaiprimex(x)\n\nScaled derivative of the Airy function of the first kind operatornameAi(x) e^frac23 x sqrtx.  Throws DomainError for negative Real arguments.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.airybi",
    "page": "List",
    "title": "SpecialFunctions.airybi",
    "category": "function",
    "text": "airybi(x)\n\nAiry function of the second kind operatornameBi(x).\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.airybiprime",
    "page": "List",
    "title": "SpecialFunctions.airybiprime",
    "category": "function",
    "text": "airybiprime(x)\n\nDerivative of the Airy function of the second kind operatornameBi(x).\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.airybix",
    "page": "List",
    "title": "SpecialFunctions.airybix",
    "category": "function",
    "text": "airybix(x)\n\nScaled Airy function of the second kind operatornameBi(x) e^- left operatornameRe left( frac23 x sqrtx right) right.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.airybiprimex",
    "page": "List",
    "title": "SpecialFunctions.airybiprimex",
    "category": "function",
    "text": "airybiprimex(x)\n\nScaled derivative of the Airy function of the second kind operatornameBi(x) e^- left operatornameRe left( frac23 x sqrtx right) right.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.besselj0",
    "page": "List",
    "title": "SpecialFunctions.besselj0",
    "category": "function",
    "text": "besselj0(x)\n\nBessel function of the first kind of order 0, J_0(x).\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.besselj1",
    "page": "List",
    "title": "SpecialFunctions.besselj1",
    "category": "function",
    "text": "besselj1(x)\n\nBessel function of the first kind of order 1, J_1(x).\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.besselj",
    "page": "List",
    "title": "SpecialFunctions.besselj",
    "category": "function",
    "text": "besselj(nu, x)\n\nBessel function of the first kind of order nu, J_nu(x).\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.besseljx",
    "page": "List",
    "title": "SpecialFunctions.besseljx",
    "category": "function",
    "text": "besseljx(nu, x)\n\nScaled Bessel function of the first kind of order nu, J_nu(x) e^-  operatornameIm(x) .\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.bessely0",
    "page": "List",
    "title": "SpecialFunctions.bessely0",
    "category": "function",
    "text": "bessely0(x)\n\nBessel function of the second kind of order 0, Y_0(x).\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.bessely1",
    "page": "List",
    "title": "SpecialFunctions.bessely1",
    "category": "function",
    "text": "bessely1(x)\n\nBessel function of the second kind of order 1, Y_1(x).\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.bessely",
    "page": "List",
    "title": "SpecialFunctions.bessely",
    "category": "function",
    "text": "bessely(nu, x)\n\nBessel function of the second kind of order nu, Y_nu(x).\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.besselyx",
    "page": "List",
    "title": "SpecialFunctions.besselyx",
    "category": "function",
    "text": "besselyx(nu, x)\n\nScaled Bessel function of the second kind of order nu, Y_nu(x) e^-  operatornameIm(x) .\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.hankelh1",
    "page": "List",
    "title": "SpecialFunctions.hankelh1",
    "category": "function",
    "text": "hankelh1(nu, x)\n\nBessel function of the third kind of order nu, H^(1)_nu(x).\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.hankelh1x",
    "page": "List",
    "title": "SpecialFunctions.hankelh1x",
    "category": "function",
    "text": "hankelh1x(nu, x)\n\nScaled Bessel function of the third kind of order nu, H^(1)_nu(x) e^-x i.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.hankelh2",
    "page": "List",
    "title": "SpecialFunctions.hankelh2",
    "category": "function",
    "text": "hankelh2(nu, x)\n\nBessel function of the third kind of order nu, H^(2)_nu(x).\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.hankelh2x",
    "page": "List",
    "title": "SpecialFunctions.hankelh2x",
    "category": "function",
    "text": "hankelh2x(nu, x)\n\nScaled Bessel function of the third kind of order nu, H^(2)_nu(x) e^x i.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.besselh",
    "page": "List",
    "title": "SpecialFunctions.besselh",
    "category": "function",
    "text": "besselh(nu, [k=1,] x)\n\nBessel function of the third kind of order nu (the Hankel function). k is either 1 or 2, selecting hankelh1 or hankelh2, respectively. k defaults to 1 if it is omitted. (See also besselhx for an exponentially scaled variant.)\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.besselhx",
    "page": "List",
    "title": "SpecialFunctions.besselhx",
    "category": "function",
    "text": "besselhx(nu, [k=1,] z)\n\nCompute the scaled Hankel function exp(iz) H_ν^(k)(z), where k is 1 or 2, H_ν^(k)(z) is besselh(nu, k, z), and  is - for k=1 and + for k=2.  k defaults to 1 if it is omitted.\n\nThe reason for this function is that H_ν^(k)(z) is asymptotically proportional to exp(iz)sqrtz for large z, and so the besselh function is susceptible to overflow or underflow when z has a large imaginary part.  The besselhx function cancels this exponential factor (analytically), so it avoids these problems.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.besseli",
    "page": "List",
    "title": "SpecialFunctions.besseli",
    "category": "function",
    "text": "besseli(nu, x)\n\nModified Bessel function of the first kind of order nu, I_nu(x).\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.besselix",
    "page": "List",
    "title": "SpecialFunctions.besselix",
    "category": "function",
    "text": "besselix(nu, x)\n\nScaled modified Bessel function of the first kind of order nu, I_nu(x) e^-  operatornameRe(x) .\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.besselk",
    "page": "List",
    "title": "SpecialFunctions.besselk",
    "category": "function",
    "text": "besselk(nu, x)\n\nModified Bessel function of the second kind of order nu, K_nu(x).\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.besselkx",
    "page": "List",
    "title": "SpecialFunctions.besselkx",
    "category": "function",
    "text": "besselkx(nu, x)\n\nScaled modified Bessel function of the second kind of order nu, K_nu(x) e^x.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.ellipk",
    "page": "List",
    "title": "SpecialFunctions.ellipk",
    "category": "function",
    "text": "ellipk(m)\n\nComputes Complete Elliptic Integral of 1st kind K(m) for parameter m given by\n\noperatornameellipk(m)\n= K(m)\n= int_0^ fracpi2  frac1sqrt1 - m sin^2 theta  mathrmdtheta\nquad textfor quad m in left( -infty 1 right  \n\nExternal links: DLMF, Wikipedia.\n\nSee also: ellipe(m).\n\nArguments\n\nm: parameter m is in relation to elliptic modulus k by k^2=m and modular   angle alpha by k=sin alpha\n\nImplementation\n\nUsing piecewise approximation polynomial as given in\n\n\'Fast Computation of Complete Elliptic Integrals and Jacobian Elliptic Functions\', Fukushima, Toshio. (2014). F09-FastEI. Celest Mech Dyn Astr, DOI 10.1007/s10569-009-9228-z, https://pdfs.semanticscholar.org/8112/c1f56e833476b61fc54d41e194c962fbe647.pdf\n\nFor m0`, followed\n\nFukushima, Toshio. (2014). \'Precise, compact, and fast computation of complete elliptic integrals by piecewise minimax rational function approximation\'. Journal of Computational and Applied Mathematics. 282. DOI 10.13140/2.1.1946.6245., https://www.researchgate.net/publication/267330394\n\nAlso suggested in this paper that we should consider domain only from (-infty1.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.ellipe",
    "page": "List",
    "title": "SpecialFunctions.ellipe",
    "category": "function",
    "text": "ellipe(m)\n\nComputes Complete Elliptic Integral of 2nd kind E(m) for parameter m given by\n\noperatornameellipe(m)\n= E(m)\n= int_0^ fracpi2  sqrt1 - m sin^2 theta  mathrmdtheta\nquad textfor quad m in left( -infty 1 right  \n\nExternal links: DLMF, Wikipedia.\n\nSee also: ellipk(m).\n\nArguments\n\nm: parameter m is in relation to elliptic modulus k by k^2=m and modular   angle alpha by k=sin alpha\n\nImplementation\n\nUsing piecewise approximation polynomial as given in\n\n\'Fast Computation of Complete Elliptic Integrals and Jacobian Elliptic Functions\', Fukushima, Toshio. (2014). F09-FastEI. Celest Mech Dyn Astr, DOI 10.1007/s10569-009-9228-z, https://pdfs.semanticscholar.org/8112/c1f56e833476b61fc54d41e194c962fbe647.pdf\n\nFor m0`, followed\n\nFukushima, Toshio. (2014). \'Precise, compact, and fast computation of complete elliptic integrals by piecewise minimax rational function approximation\'. Journal of Computational and Applied Mathematics. 282. DOI 10.13140/2.1.1946.6245., https://www.researchgate.net/publication/267330394\n\nAlso suggested in this paper that we should consider domain only from (-infty1.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.eta",
    "page": "List",
    "title": "SpecialFunctions.eta",
    "category": "function",
    "text": "eta(x)\n\nDirichlet eta function eta(s) = sum^infty_n=1(-1)^n-1n^s.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.zeta",
    "page": "List",
    "title": "SpecialFunctions.zeta",
    "category": "function",
    "text": "zeta(s, z)\n\nGeneralized zeta function zeta(s z), defined by the sum sum_k=0^infty ((k+z)^2)^-s2, where any term with k+z=0 is excluded.  For Re z  0, this definition is equivalent to the Hurwitz zeta function sum_k=0^infty (k+z)^-s.   For z=1, it yields the Riemann zeta function zeta(s).\n\n\n\n\n\nzeta(s)\n\nRiemann zeta function zeta(s).\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.gamma",
    "page": "List",
    "title": "SpecialFunctions.gamma",
    "category": "function",
    "text": "gamma(z)\n\nCompute the gamma function for complex z, defined by\n\nGamma(z)\n=\nbegincases\n    n\n     textfor quad z = n+1  n = 012dots\n    \n    int_0^infty t^z-1 mathrm e^-t  mathrm dt\n     textfor quad Re(z)  0\nendcases\n\nand by analytic continuation in the whole complex plane.\n\nExternal links: DLMF, Wikipedia.\n\nSee also: loggamma(z).\n\nImplementation by\n\nFloat: C standard math library   libm.\nComplex: by exp(loggamma(z)).\nBigFloat: C library for multiple-precision floating-point MPFR\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.gamma_inc",
    "page": "List",
    "title": "SpecialFunctions.gamma_inc",
    "category": "function",
    "text": "gamma_inc(a,x,IND)\n\nDLMF: https://dlmf.nist.gov/8.2#E4 , https://dlmf.nist.gov/8.2#E5 Wikipedia: https://en.wikipedia.org/wiki/Incompletegammafunction IND –> Accuracy desired ; IND=0 means 14 significant digits accuracy , IND=1 means 6 significant digit and IND=2 means only 3 digit accuracy suffices. gammainc(a,x) or P(a,x) is the Incomplete gamma function ratio given by : ``1/\\Gamma (a) \\int{0}^{x} e^{-t}t^{a-1} dtgamma_q(ax) or Q(ax) is the Incomplete gamma function ratio given by  1 - P(ax) -1/\\Gamma (a) \\int{x}^{\\infty} e^{-t}t^{a-1} dt`` Returns a tuple (gammainc, gammaq) where gammainc + gamma_q = 1.0\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.gamma_inc_inv",
    "page": "List",
    "title": "SpecialFunctions.gamma_inc_inv",
    "category": "function",
    "text": "gamma_inc_inv(a,p,q)\n\nDLMF: https://dlmf.nist.gov/8.2#E4 , https://dlmf.nist.gov/8.2#E5 Wiki: https://en.wikipedia.org/wiki/Incompletegammafunction\n\ngammainc(a,x) or (P(a,x),Q(a,x)) is the Incomplete gamma function ratio given by : ``1/\\Gamma (a) \\int{0}^{x} e^{-t}t^{a-1} dt`gamma_inc_inv(a,p,q) inverts the gamma_inc function, by computingxgivena,p,q` in P(a,x)=p and Q(a,x)=q.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.beta_inc",
    "page": "List",
    "title": "SpecialFunctions.beta_inc",
    "category": "function",
    "text": "beta_inc(a,b,x,y)\n\nComputes Incomplete Beta Function Ratios given by:\n\nI_x(ab) = G(ab) int_0^x t^a-1(1-t)^b-1 dt\n\nand I_y(ab) = 10 - I_x(ab). given B(ab) = 1G(ab) = Gamma(a)Gamma(b)Gamma(a+b) and x+y = 1.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.loggamma",
    "page": "List",
    "title": "SpecialFunctions.loggamma",
    "category": "function",
    "text": "loggamma(x)\n\nComputes the logarithm of gamma for given x. If x is a Real, then it throws a DomainError if gamma(x) is negative.\n\nSee also logabsgamma.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.logabsgamma",
    "page": "List",
    "title": "SpecialFunctions.logabsgamma",
    "category": "function",
    "text": "logabsgamma(x)\n\nCompute the logarithm of absolute value of gamma for Real xand returns a tuple (log(abs(gamma(x))), sign(gamma(x))).\n\nSee also loggamma.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.logfactorial",
    "page": "List",
    "title": "SpecialFunctions.logfactorial",
    "category": "function",
    "text": "logfactorial(x)\n\nCompute the logarithmic factorial of a nonnegative integer x. Equivalent to loggamma of x + 1, but loggamma extends this function to non-integer x.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.beta",
    "page": "List",
    "title": "SpecialFunctions.beta",
    "category": "function",
    "text": "beta(x, y)\n\nEuler integral of the first kind operatornameB(xy) = Gamma(x)Gamma(y)Gamma(x+y).\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.logbeta",
    "page": "List",
    "title": "SpecialFunctions.logbeta",
    "category": "function",
    "text": "logbeta(x, y)\n\nNatural logarithm of the beta function log(operatornameB(xy)).\n\nSee also logabsbeta.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.logabsbeta",
    "page": "List",
    "title": "SpecialFunctions.logabsbeta",
    "category": "function",
    "text": "logabsbeta(x, y)\n\nCompute the natural logarithm of the absolute value of the beta function, returning a tuple (log(abs(beta(x,y))), sign(beta(x,y)))\n\nSee also logbeta.\n\n\n\n\n\n"
},

{
    "location": "functions_list/#SpecialFunctions.logabsbinomial",
    "page": "List",
    "title": "SpecialFunctions.logabsbinomial",
    "category": "function",
    "text": "logabsbinomial(n, k)\n\nAccurate natural logarithm of the absolute value of the binomial coefficient binomial(n, k) for large n and k near n/2.\n\nReturns a tuple (log(abs(binomial(n,k))), sign(binomial(n,k))).\n\n\n\n\n\n"
},

{
    "location": "functions_list/#Functions-1",
    "page": "List",
    "title": "Functions",
    "category": "section",
    "text": "CurrentModule = SpecialFunctionsSpecialFunctions.erf\nSpecialFunctions.erfc\nSpecialFunctions.erfcx\nSpecialFunctions.erfi\nSpecialFunctions.dawson\nSpecialFunctions.erfinv\nSpecialFunctions.erfcinv\nSpecialFunctions.sinint\nSpecialFunctions.cosint\nSpecialFunctions.digamma\nSpecialFunctions.invdigamma\nSpecialFunctions.trigamma\nSpecialFunctions.polygamma\nSpecialFunctions.airyai\nSpecialFunctions.airyaiprime\nSpecialFunctions.airyaix\nSpecialFunctions.airyaiprimex\nSpecialFunctions.airybi\nSpecialFunctions.airybiprime\nSpecialFunctions.airybix\nSpecialFunctions.airybiprimex\nSpecialFunctions.besselj0\nSpecialFunctions.besselj1\nSpecialFunctions.besselj\nSpecialFunctions.besseljx\nSpecialFunctions.bessely0\nSpecialFunctions.bessely1\nSpecialFunctions.bessely\nSpecialFunctions.besselyx\nSpecialFunctions.hankelh1\nSpecialFunctions.hankelh1x\nSpecialFunctions.hankelh2\nSpecialFunctions.hankelh2x\nSpecialFunctions.besselh\nSpecialFunctions.besselhx\nSpecialFunctions.besseli\nSpecialFunctions.besselix\nSpecialFunctions.besselk\nSpecialFunctions.besselkx\nSpecialFunctions.ellipk\nSpecialFunctions.ellipe\nSpecialFunctions.eta\nSpecialFunctions.zeta\nSpecialFunctions.gamma\nSpecialFunctions.gamma_inc\nSpecialFunctions.gamma_inc_inv\nSpecialFunctions.beta_inc\nSpecialFunctions.loggamma\nSpecialFunctions.logabsgamma\nSpecialFunctions.logfactorial\nSpecialFunctions.beta\nSpecialFunctions.logbeta\nSpecialFunctions.logabsbeta\nSpecialFunctions.logabsbinomial"
},

]}
