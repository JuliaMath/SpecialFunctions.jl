# airy
@test_approx_eq airy(1.8) 0.0470362168668458052247
@test_approx_eq airyprime(1.8) -0.0685247801186109345638
@test_approx_eq airybi(1.8) 2.595869356743906290060
@test_approx_eq airybiprime(1.8) 2.98554005084659907283
@test_throws Base.Math.AmosException airy(200im)
@test_throws Base.Math.AmosException airybi(200)
z = 1.8 + 1.0im
@test_approx_eq airyx(0, z) airy(0, z) * exp(2/3 * z * sqrt(z))
@test_approx_eq airyx(1, z) airy(1, z) * exp(2/3 * z * sqrt(z))
@test_approx_eq airyx(2, z) airy(2, z) * exp(-abs(real(2/3 * z * sqrt(z))))
@test_approx_eq airyx(3, z) airy(3, z) * exp(-abs(real(2/3 * z * sqrt(z))))

# besselh
true_h133 = 0.30906272225525164362 - 0.53854161610503161800im
@test_approx_eq besselh(3,1,3) true_h133
@test_approx_eq besselh(-3,1,3) -true_h133
@test_approx_eq besselh(3,2,3) conj(true_h133)
@test_approx_eq besselh(-3,2,3) -conj(true_h133)
@test_throws Base.Math.AmosException besselh(1,0)


# besseli
true_i33 = 0.95975362949600785698
@test_approx_eq besseli(3,3) true_i33
@test_approx_eq besseli(-3,3) true_i33
@test_approx_eq besseli(3,-3) -true_i33
@test_approx_eq besseli(-3,-3) -true_i33
@test_throws Base.Math.AmosException besseli(1,1000)

# besselj
@test besselj(0,0) == 1
for i = 1:5
    @test besselj(i,0) == 0
    @test besselj(-i,0) == 0
end

j33 = besselj(3,3.)
@test besselj(3,3) == j33
@test besselj(-3,-3) == j33
@test besselj(-3,3) == -j33
@test besselj(3,-3) == -j33

j43 = besselj(4,3.)
@test besselj(4,3) == j43
@test besselj(-4,-3) == j43
@test besselj(-4,3) == j43
@test besselj(4,-3) == j43

@test_approx_eq j33 0.30906272225525164362
@test_approx_eq j43 0.13203418392461221033
@test_throws DomainError    besselj(0.1, -0.4)
@test_approx_eq besselj(0.1, complex(-0.4)) 0.820421842809028916 + 0.266571215948350899im
@test_approx_eq besselj(3.2, 1.3+0.6im) 0.01135309305831220201 + 0.03927719044393515275im
@test_approx_eq besselj(1, 3im) 3.953370217402609396im
@test_throws Base.Math.AmosException besselj(20,1000im)

# besselk
true_k33 = 0.12217037575718356792
@test_approx_eq besselk(3,3) true_k33
@test_approx_eq besselk(-3,3) true_k33
true_k3m3 = -0.1221703757571835679 - 3.0151549516807985776im
@test_throws DomainError besselk(3,-3)
@test_approx_eq besselk(3,complex(-3)) true_k3m3
@test_approx_eq besselk(-3,complex(-3)) true_k3m3
@test_throws Base.Math.AmosException besselk(200,0.01)
# issue #6564
@test besselk(1.0,0.0) == Inf

# bessely
y33 = bessely(3,3.)
@test bessely(3,3) == y33
@test_approx_eq bessely(-3,3) -y33
@test_approx_eq y33 -0.53854161610503161800
@test_throws DomainError bessely(3,-3)
@test_approx_eq bessely(3,complex(-3)) 0.53854161610503161800 - 0.61812544451050328724im
@test_throws Base.Math.AmosException bessely(200.5,0.1)

# issue #6653
for f in (besselj,bessely,besseli,besselk,hankelh1,hankelh2)
    @test_approx_eq f(0,1) f(0,complex128(1))
    @test_approx_eq f(0,1) f(0,complex64(1))
end

# scaled bessel[ijky] and hankelh[12]
for x in (1.0, 0.0, -1.0), y in (1.0, 0.0, -1.0), nu in (1.0, 0.0, -1.0)
    z = complex128(x + y * im)
    z == zero(z) || @test_approx_eq hankelh1x(nu, z) hankelh1(nu, z) * exp(-z * im)
    z == zero(z) || @test_approx_eq hankelh2x(nu, z) hankelh2(nu, z) * exp(z * im)
    (nu < 0 && z == zero(z)) || @test_approx_eq besselix(nu, z) besseli(nu, z) * exp(-abs(real(z)))
    (nu < 0 && z == zero(z)) || @test_approx_eq besseljx(nu, z) besselj(nu, z) * exp(-abs(imag(z)))
    z == zero(z) || @test_approx_eq besselkx(nu, z) besselk(nu, z) * exp(z)
    z == zero(z) || @test_approx_eq besselyx(nu, z) bessely(nu, z) * exp(-abs(imag(z)))
end
@test_throws Base.Math.AmosException hankelh1x(1, 0)
@test_throws Base.Math.AmosException hankelh2x(1, 0)
@test_throws Base.Math.AmosException besselix(-1, 0)
@test_throws Base.Math.AmosException besseljx(-1, 0)
@test besselkx(1, 0) == Inf
@test_throws Base.Math.AmosException besselyx(1, 0)
