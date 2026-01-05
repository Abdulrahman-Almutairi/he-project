#pragma once
#include "bfv/poly.h"

namespace bfv
{
  struct SecretKey
  {
    Poly s;
  };

  struct PublicKey
  {
    Poly b;
    Poly a;
  };

  struct KeyPair
  {
    PublicKey pk;
    SecretKey sk;
  };
}
