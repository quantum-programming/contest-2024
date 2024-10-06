import numbers
from functools import cached_property, total_ordering
import mpmath

from .mymath import SQRT2, sign, floorsqrt, rounddiv, pow_sqrt2


@total_ordering
class ZRootTwo():
    def __init__(self, a, b):
        self._a = a
        self._b = b

    @property
    def a(self):
        return self._a

    @property
    def b(self):
        return self._b

    @cached_property
    def coef(self):
        return [self._a, self._b]

    def __repr__(self):
        return f"ZRootTwo({self._a}, {self._b})"

    def __str__(self):
        return f"{self._a}{self._b:+}√2"

    @classmethod
    def fromInt(cls, x):
        return cls(x, 0)

    @classmethod
    def fromZOmega(cls, x):
        if x.b == 0 and x.a == -x.c:
            return cls(x.d, x.c)
        else:
            raise ValueError

    def __eq__(self, other):
        if isinstance(other, numbers.Integral):
            return self == self.fromInt(other)
        elif isinstance(other, self.__class__):
            return self.a == other.a and self.b == other.b
        else:
            return False

    def __lt__(self, other):
        if isinstance(other, numbers.Integral):
            return self < self.fromInt(other)
        elif isinstance(other, self.__class__):
            if self.b < other.b:
                return self.a < other.a or (self.a - other.a) ** 2 < 2 * (self.b - other.b) ** 2
            else:
                return self.a < other.a and (self.a - other.a) ** 2 > 2 * (self.b - other.b) ** 2
        else:
            return False

    def __add__(self, other):
        if isinstance(other, numbers.Integral):
            return self + self.fromInt(other)
        elif isinstance(other, self.__class__):
            return self.__class__(self._a + other.a, self._b + other.b)
        else:
            return NotImplemented

    def __radd__(self, other):
        if isinstance(other, numbers.Integral):
            return self + other
        elif isinstance(other, self.__class__):
            return self + other
        else:
            return NotImplemented

    def __sub__(self, other):
        return self + (- other)

    def __rsub__(self, other):
        return (-self) + other

    def __neg__(self):
        return self.__class__(-self._a, -self._b)

    def __mul__(self, other):
        if isinstance(other, numbers.Integral):
            return self * self.fromInt(other)
        elif isinstance(other, self.__class__):
            new_a = self._a * other.a + 2 * self._b * other.b
            new_b = self._a * other.b + self._b * other.a
            return self.__class__(new_a, new_b)
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, numbers.Integral):
            return self * other
        elif isinstance(other, self.__class__):
            return self * other
        else:
            return NotImplemented

    @cached_property
    def inv(self):
        if self.norm == 1:
            return self.conj_sq2
        elif self.norm == -1:
            return -self.conj_sq2
        else:
            raise ZeroDivisionError

    def __pow__(self, other):
        if isinstance(other, numbers.Integral):
            if other < 0:
                return self.inv ** -other
            else:
                new = self.__class__(1, 0)
                tmp = self
                while other > 0:
                    if other & 1:
                        new *= tmp
                    tmp *= tmp
                    other >>= 1
                return new
        else:
            return NotImplemented

    def sqrt(self):
        norm = self.norm
        if norm < 0 or self._a < 0:
            return None
        r = floorsqrt(norm)
        a1 = floorsqrt((self._a + r) // 2)
        b1 = floorsqrt((self._a - r) // 4)
        a2 = floorsqrt((self._a - r) // 2)
        b2 = floorsqrt((self._a + r) // 4)
        if sign(self._a) * sign(self._b) >= 0:
            w1 = ZRootTwo(a1, b1)
            w2 = ZRootTwo(a2, b2)
        else:
            w1 = ZRootTwo(a1, -b1)
            w2 = ZRootTwo(a2, -b2)
        if self == w1 * w1:
            return w1
        elif self == w2 * w2:
            return w2
        else:
            return None

    def __divmod__(self, other):
        if isinstance(other, numbers.Integral):
            return divmod(self, self.fromInt(other))
        elif isinstance(other, self.__class__):
            p = self * other.conj_sq2
            k = other.norm
            q = self.__class__(rounddiv(p.a, k), rounddiv(p.b, k))
            r = self - other * q
            return q, r
        else:
            return NotImplemented

    def __rdivmod__(self, other):
        if isinstance(other, numbers.Integral):
            return divmod(self.fromInt(other), self)
        elif isinstance(other, self.__class__):
            return divmod(other, self)
        else:
            return NotImplemented

    def __floordiv__(self, other):
        q, _ = divmod(self, other)
        return q

    def __rfloordiv__(self, other):
        q, _ = divmod(other, self)
        return q

    def __mod__(self, other):
        _, r = divmod(self, other)
        return r

    def __rmod__(self, other):
        _, r = divmod(other, self)
        return r

    @classmethod
    def sim(cls, a, b):
        return a % b == 0 and b % a == 0

    @classmethod
    def ext_gcd(cls, a, b):
        if isinstance(a, numbers.Integral):
            a = cls.fromInt(a)
        if isinstance(b, numbers.Integral):
            b = cls.fromInt(b)
        x = cls.fromInt(1)
        y = cls.fromInt(0)
        z = cls.fromInt(0)
        w = cls.fromInt(1)
        while b != 0:
            q, r = divmod(a, b)
            x, y = y, x - y * q
            z, w = w, z - w * q
            a, b = b, r
        return x, z, a

    @classmethod
    def gcd(cls, a, b):
        _, _, g = cls.ext_gcd(a, b)
        return g

    @cached_property
    def parity(self):
        return self._a & 1

    @cached_property
    def norm(self):
        return self._a ** 2 - 2 * self._b ** 2

    @cached_property
    def toReal(self):
        return self._a + SQRT2() * self._b

    @cached_property
    def conj_sq2(self):
        return self.__class__(self._a, -self._b)


@total_ordering
class DRootTwo():
    def __init__(self, alpha, k):
        self._alpha = alpha
        self._k = k

    @property
    def alpha(self):
        return self._alpha

    @property
    def k(self):
        return self._k

    def __repr__(self):
        return f"DRootTwo({self._alpha}, {self._k})"

    def __str__(self):
        return f"{self._alpha} / √2^{self._k}"

    @classmethod
    def fromInt(cls, x):
        return cls(ZRootTwo.fromInt(x), 0)

    @classmethod
    def fromZRootTwo(cls, x):
        return cls(x, 0)

    @classmethod
    def fromZOmega(cls, x):
        return cls(ZRootTwo.fromZOmega(x), 0)

    @classmethod
    def fromDOmega(cls, x):
        return cls(ZRootTwo.fromZOmega(x.u), x.k)

    def renew_denomexp(self, new_k):
        new_alpha = self.mul_sqrt2(new_k - self.k).alpha
        return self.__class__(new_alpha, new_k)

    def mul_sqrt2(self, d):
        if d < 0:
            d_div_2, d_mod_2 = (-d) >> 1, (-d) & 1
            if d_mod_2 == 0:
                bit = (1 << d_div_2) - 1
                if self._alpha.a & bit == 0 and self._alpha.b & bit == 0:
                    new_alpha = ZRootTwo(self._alpha.a >> d_div_2, self._alpha.b >> d_div_2)
                else:
                    raise ValueError
            else:
                bit = (1 << d_div_2) - 1
                bit2 = (1 << (d_div_2 + 1)) - 1
                if self._alpha.a & bit2 == 0 and self._alpha.b & bit == 0:
                    new_alpha = ZRootTwo(self._alpha.b >> d_div_2, self._alpha.a >> (d_div_2 + 1))
                else:
                    raise ValueError
            return self.__class__(new_alpha, self.k)
        else:
            d_div_2, d_mod_2 = d >> 1, d & 1
            new_alpha = self._alpha * (1 << d_div_2)
            if d_mod_2:
                new_alpha *= ZRootTwo(0, 1)
            return self.__class__(new_alpha, self.k)

    def mul_sqrt2_renewing_denomexp(self, d):
        if d > self._k:
            raise ValueError
        return self.__class__(self._alpha, self._k - d)

    def __eq__(self, other):
        if isinstance(other, numbers.Integral):
            return self == self.fromInt(other)
        elif isinstance(other, ZRootTwo):
            return self == self.fromZRootTwo(other)
        elif isinstance(other, self.__class__):
            if self._k < other.k:
                return self.renew_denomexp(other.k) == other
            elif self._k > other.k:
                return self == other.renew_denomexp(self._k)
            else:
                return self._alpha == other.alpha and self._k == other.k
        else:
            return False

    def __lt__(self, other):
        if isinstance(other, numbers.Integral):
            return self < self.fromInt(other)
        elif isinstance(other, ZRootTwo):
            return self < self.fromZRootTwo(other)
        elif isinstance(other, self.__class__):
            if self._k < other.k:
                return self.renew_denomexp(other.k) < other
            elif self._k > other.k:
                return self < other.renew_denomexp(self._k)
            else:
                return self._alpha < other.alpha
        else:
            return False

    def __add__(self, other):
        if isinstance(other, numbers.Integral):
            return self + self.fromInt(other)
        elif isinstance(other, ZRootTwo):
            return self + self.fromZRootTwo(other)
        elif isinstance(other, self.__class__):
            if self._k < other.k:
                return self.renew_denomexp(other.k) + other
            elif self._k > other.k:
                return self + other.renew_denomexp(self._k)
            else:
                return self.__class__(self._alpha + other.alpha, self._k)
        else:
            return NotImplemented

    def __radd__(self, other):
        if isinstance(other, numbers.Integral):
            return self + other
        elif isinstance(other, ZRootTwo):
            return self + other
        elif isinstance(other, self.__class__):
            return self + other
        else:
            return NotImplemented

    def __sub__(self, other):
        return self + (- other)

    def __rsub__(self, other):
        return (-self) + other

    def __neg__(self):
        return self.__class__(-self._alpha, self._k)

    def __mul__(self, other):
        if isinstance(other, numbers.Integral):
            return self * self.fromInt(other)
        elif isinstance(other, ZRootTwo):
            return self * self.fromZRootTwo(other)
        elif isinstance(other, self.__class__):
            return self.__class__(self._alpha * other.alpha, self._k + other.k)
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, numbers.Integral):
            return self * other
        elif isinstance(other, ZRootTwo):
            return self * other
        elif isinstance(other, self.__class__):
            return self * other
        else:
            return NotImplemented

    @cached_property
    def parity(self):
        return self._alpha.parity

    @cached_property
    def scale(self):
        return pow_sqrt2(self._k)

    @cached_property
    def scale_pow_2(self):
        return 1 << self._k

    @cached_property
    def toReal(self):
        return self.alpha.toReal / self.scale

    @cached_property
    def conj_sq2(self):
        return (self.__class__(- self._alpha.conj_sq2, self._k) if self._k & 1
                else self.__class__(self._alpha.conj_sq2, self._k))

    @classmethod
    def pow_sqrt2_inv(cls, k):
        return cls(ZRootTwo(1, 0), k)


class ZOmega():
    def __init__(self, a, b, c, d):
        self._a = a
        self._b = b
        self._c = c
        self._d = d

    @property
    def a(self):
        return self._a

    @property
    def b(self):
        return self._b

    @property
    def c(self):
        return self._c

    @property
    def d(self):
        return self._d

    @cached_property
    def coef(self):
        return [self._d, self._c, self._b, self._a]

    def __repr__(self):
        return f"ZOmega({self._a}, {self._b}, {self._c}, {self._d})"

    def __str__(self):
        return f"{self._a}ω^3{self._b:+}ω^2{self._c:+}ω{self._d:+}"

    @classmethod
    def fromInt(cls, x):
        return cls(0, 0, 0, x)

    @classmethod
    def fromZRootTwo(cls, x):
        return cls(-x.b, 0, x.b, x.a)

    def __eq__(self, other):
        if isinstance(other, numbers.Integral):
            return self == self.fromInt(other)
        elif isinstance(other, ZRootTwo):
            return self == self.fromZRootTwo(other)
        elif isinstance(other, self.__class__):
            return (self.a == other.a and self.b == other.b
                    and self.c == other.c and self.d == other.d)
        else:
            return False

    def __add__(self, other):
        if isinstance(other, numbers.Integral):
            return self + self.fromInt(other)
        elif isinstance(other, ZRootTwo):
            return self + self.fromZRootTwo(other)
        elif isinstance(other, self.__class__):
            return self.__class__(self._a + other.a, self._b + other.b,
                                  self._c + other.c, self._d + other.d)
        else:
            return NotImplemented

    def __radd__(self, other):
        if isinstance(other, numbers.Integral):
            return self + other
        elif isinstance(other, ZRootTwo):
            return self + other
        elif isinstance(other, self.__class__):
            return self + other
        else:
            return NotImplemented

    def __sub__(self, other):
        return self + (- other)

    def __rsub__(self, other):
        return (-self) + other

    def __neg__(self):
        return self.__class__(-self._a, -self._b, -self._c, -self._d)

    def __mul__(self, other):
        if isinstance(other, numbers.Integral):
            return self * self.fromInt(other)
        elif isinstance(other, ZRootTwo):
            return self * self.fromZRootTwo(other)
        elif isinstance(other, self.__class__):
            new_coef = [0] * 4
            for i in range(4):
                for j in range(4):
                    if i + j < 4:
                        new_coef[i + j] += self.coef[i] * other.coef[j]
                    else:
                        new_coef[i + j - 4] -= self.coef[i] * other.coef[j]
            return self.__class__(*reversed(new_coef))
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, numbers.Integral):
            return self * other
        elif isinstance(other, ZRootTwo):
            return self * other
        elif isinstance(other, self.__class__):
            return self * other
        else:
            return NotImplemented

    @cached_property
    def inv(self):
        if self.norm == 1:
            return self.conj_sq2 * self.conj * self.conj.conj_sq2
        else:
            raise ZeroDivisionError

    def __pow__(self, other):
        if isinstance(other, int):
            if other < 0:
                return NotImplemented
            else:
                new = self.fromInt(1)
                tmp = self
                while other > 0:
                    if other & 1:
                        new *= tmp
                    tmp *= tmp
                    other >>= 1
                return new
        else:
            return NotImplemented

    def __divmod__(self, other):
        if isinstance(other, numbers.Integral):
            return divmod(self, self.fromInt(other))
        elif isinstance(other, ZRootTwo):
            return divmod(self, self.fromZRootTwo(other))
        elif isinstance(other, self.__class__):
            p = self * other.conj * other.conj.conj_sq2 * other.conj_sq2
            k = other.norm
            q = self.__class__(rounddiv(p.a, k), rounddiv(p.b, k),
                               rounddiv(p.c, k), rounddiv(p.d, k))
            r = self - other * q
            return q, r
        else:
            return NotImplemented

    def __rdivmod__(self, other):
        if isinstance(other, numbers.Integral):
            return divmod(self.fromInt(other), self)
        elif isinstance(other, ZRootTwo):
            return divmod(self.fromZRootTwo(other), self)
        elif isinstance(other, self.__class__):
            return divmod(other, self)
        else:
            return NotImplemented

    def __floordiv__(self, other):
        q, _ = divmod(self, other)
        return q

    def __rfloordiv__(self, other):
        q, _ = divmod(other, self)
        return q

    def __mod__(self, other):
        _, r = divmod(self, other)
        return r

    def __rmod__(self, other):
        _, r = divmod(other, self)
        return r

    @classmethod
    def sim(cls, a, b):
        return a % b == 0 and b % a == 0

    @classmethod
    def ext_gcd(cls, a, b):
        if isinstance(a, numbers.Integral):
            a = cls.fromInt(a)
        elif isinstance(a, ZRootTwo):
            a = cls.fromZRootTwo(a)
        if isinstance(b, numbers.Integral):
            b = cls.fromInt(b)
        elif isinstance(b, ZRootTwo):
            b = cls.fromZRootTwo(b)
        x = cls.fromInt(1)
        y = cls.fromInt(0)
        z = cls.fromInt(0)
        w = cls.fromInt(1)
        while b != 0:
            q, r = divmod(a, b)
            x, y = y, x - y * q
            z, w = w, z - w * q
            a, b = b, r
        return x, z, a

    @classmethod
    def gcd(cls, a, b):
        _, _, g = cls.ext_gcd(a, b)
        return g

    @cached_property
    def norm(self):
        return ((self._a ** 2 + self._b ** 2 + self._c ** 2 + self._d ** 2) ** 2
                - 2 * (self._a * self._b + self._b * self._c
                       + self._c * self._d - self._d * self._a) ** 2)

    @cached_property
    def real(self):
        return self._d + SQRT2() * (self._c - self._a) / 2

    @cached_property
    def imag(self):
        return self._b + SQRT2() * (self._c + self._a) / 2

    @cached_property
    def toComplex(self):
        return self.real + 1.j * self.imag

    @cached_property
    def toVector(self):
        return mpmath.matrix([self.real, self.imag])

    @cached_property
    def conj(self):
        return self.__class__(-self._c, -self._b, -self._a, self._d)

    @cached_property
    def conj_sq2(self):
        return self.__class__(-self._a, self._b, -self._c, self._d)


class DOmega():
    def __init__(self, u, k):
        self._u = u
        self._k = k

    @property
    def u(self):
        return self._u

    @property
    def k(self):
        return self._k

    def __repr__(self):
        return f"DOmega({repr(self._u)}, {self._k})"

    def __str__(self):
        return f"{self._u} / √2^{self._k}"

    @classmethod
    def fromInt(cls, x):
        return cls(ZOmega.fromInt(x), 0)

    @classmethod
    def fromZRootTwo(cls, x):
        return cls(ZOmega.fromZRootTwo(x), 0)

    @classmethod
    def fromDRootTwo(cls, x):
        return cls(ZOmega.fromZRootTwo(x.alpha), x.k)

    @classmethod
    def fromDRootTwoVector(cls, x, y, k):
        return (cls.fromDRootTwo(x) + cls.fromDRootTwo(y) * ZOmega(0, 1, 0, 0)).renew_denomexp(k)

    @classmethod
    def fromZOmega(cls, x):
        return cls(x, 0)

    def renew_denomexp(self, new_k):
        new_u = self.mul_sqrt2(new_k - self.k).u
        return self.__class__(new_u, new_k)

    def mul_sqrt2(self, d):
        if d < 0:
            d_div_2, d_mod_2 = (-d) >> 1, (-d) & 1
            if d_mod_2 == 0:
                bit = (1 << d_div_2) - 1
                if (self._u.a & bit == 0 and self._u.b & bit == 0
                        and self._u.c & bit == 0 and self._u.c & bit == 0):
                    new_u = ZOmega(self._u.a >> d_div_2, self._u.b >> d_div_2,
                                   self._u.c >> d_div_2, self._u.d >> d_div_2)
                else:
                    raise ValueError
            else:
                bit = (1 << (d_div_2 + 1)) - 1
                if ((self._u.b - self._u.d) & bit == 0
                        and (self._u.c + self._u.a) & bit == 0
                        and (self._u.b + self._u.d) & bit == 0
                        and (self._u.c - self._u.a) & bit == 0):
                    new_u = ZOmega((self._u.b - self._u.d) >> (d_div_2 + 1),
                                   (self._u.c + self._u.a) >> (d_div_2 + 1),
                                   (self._u.b + self._u.d) >> (d_div_2 + 1),
                                   (self._u.c - self._u.a) >> (d_div_2 + 1))
                else:
                    raise ValueError
            return self.__class__(new_u, self.k)
        else:
            d_div_2, d_mod_2 = d >> 1, d & 1
            new_u = self._u * (1 << d_div_2)
            if d_mod_2:
                new_u *= ZOmega(-1, 0, 1, 0)
            return self.__class__(new_u, self.k)

    def __eq__(self, other):
        if isinstance(other, numbers.Integral):
            return self == self.fromInt(other)
        elif isinstance(other, ZRootTwo):
            return self == self.fromZRootTwo(other)
        elif isinstance(other, ZOmega):
            return self == self.fromZOmega(other)
        elif isinstance(other, self.__class__):
            if self._k < other.k:
                return self.renew_denomexp(other.k) == other
            elif self._k > other.k:
                return self == other.renew_denomexp(self._k)
            else:
                return self._u == other.u and self._k == other.k
        else:
            return False

    def __add__(self, other):
        if isinstance(other, numbers.Integral):
            return self + self.fromInt(other)
        elif isinstance(other, ZRootTwo):
            return self + self.fromZRootTwo(other)
        elif isinstance(other, ZOmega):
            return self + self.fromZOmega(other)
        elif isinstance(other, self.__class__):
            if self._k < other.k:
                return self.renew_denomexp(other.k) + other
            elif self._k > other.k:
                return self + other.renew_denomexp(self._k)
            else:
                return self.__class__(self._u + other.u, self._k)
        else:
            return NotImplemented

    def __radd__(self, other):
        if isinstance(other, numbers.Integral):
            return self + other
        elif isinstance(other, ZRootTwo):
            return self + other
        elif isinstance(other, ZOmega):
            return self + other
        elif isinstance(other, self.__class__):
            return self + other
        else:
            return NotImplemented

    def __sub__(self, other):
        return self + (- other)

    def __rsub__(self, other):
        return (-self) + other

    def __neg__(self):
        return self.__class__(-self._u, self._k)

    def __mul__(self, other):
        if isinstance(other, numbers.Integral):
            return self * self.fromInt(other)
        elif isinstance(other, ZRootTwo):
            return self * self.fromZRootTwo(other)
        elif isinstance(other, ZOmega):
            return self * self.fromZOmega(other)
        elif isinstance(other, self.__class__):
            return self.__class__(self._u * other.u, self._k + other.k)
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, numbers.Integral):
            return self * other
        elif isinstance(other, ZRootTwo):
            return self * other
        elif isinstance(other, ZOmega):
            return self * other
        elif isinstance(other, self.__class__):
            return self * other
        else:
            return NotImplemented

    @cached_property
    def scale(self):
        return pow_sqrt2(self._k)

    @cached_property
    def scale_pow_2(self):
        return 1 << self._k

    @cached_property
    def real(self):
        return self._u.real / self.scale

    @cached_property
    def imag(self):
        return self._u.imag / self.scale

    @cached_property
    def toComplex(self):
        return self.real + 1.j * self.imag

    @cached_property
    def toVector(self):
        return mpmath.matrix([self.real, self.imag])

    @cached_property
    def conj(self):
        return self.__class__(self._u.conj, self._k)

    @cached_property
    def conj_sq2(self):
        return (self.__class__(- self._u.conj_sq2, self._k) if self._k & 1
                else self.__class__(self._u.conj_sq2, self._k))


LAMBDA = ZRootTwo(1, 1)
