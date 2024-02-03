# Which is faster?

from timeit import timeit

number = 10000
ideas = (
    'where(n<0, nan, sqrt(abs(n)))',
    'where(n<0, nan, sqrt(where(n<0, 0, n)))',
    'n ** where(n<0, nan, 0.5)',
)

print('Scalar:')  # Single where() is faster.

for idea in ideas:
    print(timeit(
        idea,
        number=number,
        setup=(
            'from numpy import nan, sqrt, where;'
            'n = -5.0'
        ),
    ))

print('Big array:')  # Where + abs is fastest.

for idea in ideas:
    print(timeit(
        idea,
        number=number,
        setup=(
            'from numpy import arange, nan, sqrt, where;'
            'n = - arange(5000)'
        ),
    ))
