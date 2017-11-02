from random import randint
from math import sqrt
from functools import reduce
from numpy import linalg, array

class Point(object):
    
    min_x1, max_x1 = -15, 30
    min_x2, max_x2 = -35, 15
    x10, x20 = (min_x1 + max_x1) / 2, (min_x2 + max_x2) / 2
    delta_x1, delta_x2 = abs(min_x1 - max_x1) / 2, abs(min_x2 - max_x2) / 2


    def __init__(self, n, x1, x2):
        self.n = n
        self.x1 = x1
        self.x2 = x2

        self.responses = []
        self._dispersion = 0

    def update(self):
        self._average = sum(self.responses) / len(self.responses)
        self._dispersion = sum([(p - self._average)**2 for p in self.responses]) / len(self.responses)
    
    @property
    def dispersion(self):
        return self._dispersion

    @property
    def average(self):
        return self._average

    def naturalize(self):
        self.x1 = Point.x10 + Point.delta_x1 * self.x1
        self.x2 = Point.x20 + Point.delta_x2 * self.x2

    def row(self):
        return '{0}\t{1}\t{2}\t{3}'.format(self.n, self.x1, self.x2, self.average)


def romanovsky_test():
    global plan, PROBABILITY
    
    responses_count = len(plan[0].responses)
    print('\t{} trials:'.format(responses_count))
    romanovsky_test_CRITICAL = {
        0.99: [1.73, 1.8375, 1.945, 2.0525, 2.16, 2.295, 2.43, 2.5250000000000004, 2.62, 2.685, 2.75, 2.8, 2.85, 2.9, 2.936, 2.972, 3.008, 3.044, 3.08],
        0.98: [1.72, 1.8225, 1.9249999999999998, 2.0275, 2.13, 2.25, 2.37, 2.455, 2.54, 2.6, 2.66, 2.7066666666666666, 2.7533333333333334, 2.8, 2.832, 2.864, 2.896, 2.928, 2.96],
        0.95: [1.71, 1.8075, 1.905, 2.0025, 2.1, 2.185, 2.27, 2.34, 2.41, 2.465, 2.52, 2.56, 2.6, 2.64, 2.668, 2.696, 2.7239999999999998, 2.752, 2.78],
        0.90: [1.69, 1.7675, 1.845, 1.9224999999999999, 2.0, 2.085, 2.17, 2.23, 2.29, 2.34, 2.39, 2.4233333333333333, 2.456666666666667, 2.49, 2.516, 2.5420000000000003, 2.568, 2.5940000000000003, 2.62]
    }[PROBABILITY][responses_count-2]

    fundeviation = sqrt(abs(2 * (2 * responses_count - 2) / responses_count / (responses_count - 4)))
    for i in range(len(plan)):
        for j in range(len(plan)):
            F_uv = plan[i].dispersion / plan[j].dispersion 
            if F_uv < 1:
                F_uv = 1 / F_uv
            if abs((responses_count - 2) * F_uv / responses_count - 1) / fundeviation >= romanovsky_test_CRITICAL:
                print('\t\tRomanovsky test failed')
                return False

    print('\t\tRomanovsky test passed\n')
    return True

plan = [
    Point(1, -1, -1),
    Point(2, -1, +1),
    Point(3, +1, -1)
]

variant = 13
PROBABILITY = 0.99

response_range = ((20 - variant) * 10,
                  (30 - variant) * 10)

for point in plan:
    point.responses.extend([randint(*response_range) for _ in range(1)])

print('Conducting an experiment\n')
while(True):
    for point in plan:
        point.responses.append(randint(*response_range))
        point.update()
    print('Checking HOV:')
    if romanovsky_test():
        break


m_x1, m_x2, a1, a2, a3, a11, a22, m_y = \
    [sum(x) / len(plan) for x in zip(*map(lambda p: (p.x1, 
                                                     p.x2, 
                                                     p.x1**2, 
                                                     p.x1*p.x2, 
                                                     p.x2**2, 
                                                     p.x1*p.average, 
                                                     p.x2*p.average, 
                                                     p.average
                                                                    ), plan))]

sle = [[   1, m_x1, m_x2],
       [m_x1,   a1,   a2],
       [m_x2,   a2,   a3]]


b0, b1, b2 = linalg.solve(array(sle), array([m_y, a11, a22]))



a0 = b0 - b1 * Point.x10 / Point.delta_x1 - b2 * Point.x20 / Point.delta_x2
a1 = b1 / Point.delta_x1
a2 = b2 / Point.delta_x2


print('N\tx1\t\tx2\t\tm(y)\t\ttheoretical responses')
for p in plan:
    p.naturalize()
    print(p.row() + '\t\t{0:.2f}'.format(a0 + a1 * p.x1 + a2 * p.x2))

print('\nRegression results are correct\n')

equation = '{3} equation: y = {0:.2f} + {1:.2f}*x1 + {2:.2f}*x2\n'

print(equation.format(b0, b1, b2, 'Normalized'))

print(equation.format(a0, a1, a2, 'Naturalized'))