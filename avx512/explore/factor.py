from sage.all import *
from gmpy2 import powmod, gcd
import random
from typing import Any, Callable
from sage.all import *
import random
from time import time
from gmpy2 import mpz, powmod # type: ignore

# Speed up key generation
proof.arithmetic(False)

precomputed_primes = {
    1536: (
        2294317681846666929032827582915871147244238987997154260752547927761680407421350229668773012083293163939930095650000966848974846426935880157396423658936168303309628291888369175257564190205921538136908065555094603890577067744062235964934772618271361628841709528345147510719656157699030950075087622569709218424737738471399638675144865007958387674395112624472783726720735207890567612000139879052278793348059850324662301473079411528931593319130299609354852366107441239,
        1802494200593124212877152220232057969820149101914917780567591041508929337113672371831686701300162645584692094587236891084394805966156657242794183368629797943879758663215081196373567865079353053014764237157020467354213814092738224751701006694341238296353856053396495458533514247995304409039019527215364796433416900272253602642367831373391296602437921809825075367037815161801942356870512197084830383932961071437649554131171040771475751447685105744841506905509803343
    ),
    512: (
        8928904469020647370148152532849395837667323616592943921195370500300314783770275217852294994305845832186344745870582852663110324209358748993742361744147667,
        10656148163168140310303511988105075679233093468014936755632012488407409042342935618188287095648180864781921283254707199748626066826006485182891632456057703
    ),
    768: (
        1329556492381818243302547678286484152560557709901688402333693896063703498923056705804721832732849414283887027408445712587233557135984663043747304347541113331802449222095801963080060653306593757732467353409282511419836980781555892999,
        1001398871047444817351121606809414503432013557991696786622593224176825363990201139658879587272219016290334996128275535452945154007518276707036614960005388284491331144963791645262179970088825744680279855279893302416648578043982438307
    ),
    1024: (
        172306198491877239153462054098637752223712082287916652463692759507763026075543785182048534359253307268441645596067250536115474155273834308725831759445808654026821098006583801711524840823316410138931863782308570866828886597744911362245477752845773079762297487169791854569323933861546273839610451585118589147727,
        141900624817711679032437826229884277276351993491704219282089648619873922155269135622838861926407160074876164333864883774747658092839069136111110085925902889961918040993715662225750763646988899503828356003794665525262002182119065945155914547325786676823347251109574193636418824061372987906299620546537199867527
    )
}

def generate_safe_prime(bits):
    # Very slow:
    p = 9
    while not is_prime((p - 1) // 2):
        p = random_prime(2**bits, lbound=2**(bits - 1))
    return p

def keygen(prime_len=1536):
    """Outputs a public key (n, g) and a private key (lamb, mu) for the Paillier cryptosystem."""
    if prime_len in precomputed_primes:
        p, q = precomputed_primes[prime_len]
    else:
        p = generate_safe_prime(prime_len)
        q = p
        while p == q:
            q = generate_safe_prime(prime_len)
    n = p * q
    return (mpz(p), mpz(q)), mpz(n)

def generate_generator_sagemath(n, w, g_0: int | None = None):
    """Generates a generator of order phi(N) / 4 with high probability."""
    if g_0 is None:
        g_0 = randrange(0, n ** (w + 1))
    # SageMath is dog slow:
    g = power_mod(g_0, 2 * n ** w, n ** (w + 1))
    # Equally slow:
    # g = (g_0 * g_0) % n ** (w + 1)
    # for _ in range(w):
    #     g = power_mod(g, n, n ** (w + 1))
    return g

def generate_generator(n, w, g_0: int | None = None):
    """Generates a generator of order phi(N) / 4 with high probability."""
    n_to_w = n ** w
    n_to_w_1 = n_to_w * n
    if g_0 is None:
        g_0 = random.randrange(0, n_to_w_1)
    g = powmod(g_0, 2 * n_to_w, n_to_w_1)
    return g
    # SageMath is dog slow:
    # g = power_mod(g_0, 2 * n_to_w, n_to_w_1)
    # Equally slow:
    # g = (g_0 * g_0) % n_to_w_1
    # for _ in range(w):
    #     g = power_mod(g, n, n_to_w_1)

def generate_generator_crt(p, q, n, w, g_0: int | None = None):
    """Generates a generator of order phi(N) / 4 with high probability."""
    # TODO: use CRT optimization to reduce modulus mod p and q

def benchmark_gen_generator(generate_generator_func: Callable[[int, int], Any], prime_len: int, w: int):
    print(f'Benchmarking {generate_generator_func.__name__}. {prime_len = }, {w = }')
    (p, q), n = keygen(prime_len)
    print('Keygen complete')
    start = time()
    generate_generator_func(n, w)
    end = time()
    print(f'{generate_generator_func.__name__}({prime_len = }, {w = }) Time taken:', end - start)

def test_gen_generator(generate_generator_func: Callable[[int, int, int], Any], prime_len: int, w: int):
    print(f'Testing {generate_generator_func.__name__}. {prime_len = }, {w = }')
    (p, q), n = keygen(prime_len)
    g_0 = randrange(0, n ** (w + 1))
    g = generate_generator_func(n, w, g_0)
    g_sagemath = generate_generator_sagemath(n, w, g_0)
    assert g == g_sagemath, f'Expected {g_sagemath}, got {g}'
    print('Test passed')

#### ^ copied from common.py

(p, q), N = keygen(512)
print('p =', p)
print('q =', q)
print('N =', N)

r = random_prime(2**512, lbound=2**511)
challenge = r * (p - 1) * (q - 1) // 4

for _ in range(10):
    w = random.randrange(N)
    res = powmod(w, challenge, N)
    if res == 1 or res == N - 1:
        print('Fail')
    else:
        print('Success')
        # print(res)
        found_p = gcd(res - 1, N)
        found_q = N // found_p
        print(found_p)
        print(found_q)
