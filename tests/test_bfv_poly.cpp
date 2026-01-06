#include "bfv/params.h"
#include "bfv/poly.h"
#include <cassert>
#include <iostream>

static void test_poly_add_sub_mul()
{
    bfv::Params P{1024, 12289, 2}; // small demo: N=8, q prime=97, t=17
    bfv::Poly a(P), b(P);

    a[0] = 3;
    a[1] = 5;
    a[2] = 7;
    b[0] = 10;
    b[1] = 20;
    b[2] = 30;

    auto c = bfv::add(a, b);
    assert(c[0] == (3 + 10) % 97);
    assert(c[1] == (5 + 20) % 97);
    assert(c[2] == (7 + 30) % 97);

    auto d = bfv::sub(c, b);
    assert(d[0] == a[0]);
    assert(d[1] == a[1]);
    assert(d[2] == a[2]);

    // multiply by 1 polynomial
    bfv::Poly one(P);
    one[0] = 1;
    auto e = bfv::mul_negacyclic_schoolbook(a, one);
    for (std::size_t i = 0; i < P.N; ++i)
        assert(e[i] == a[i]);

    std::cout << "BFV Poly add/sub/mul (schoolbook) passed.\n";
}

int main()
{
    test_poly_add_sub_mul();
    return 0;
}
