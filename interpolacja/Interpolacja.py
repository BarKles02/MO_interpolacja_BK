from gauss import *


class Interpolation:

    def __init__(self, points, derX, derY):
        self.points = points
        self.diffQuotients = self.computeDifferenceQuotients()
        self.h = points[1][0] - points[0][0]
        self.progressiveDiffs = self.computeProgressiveDifferences()
        self.derX = derX
        self.derY = derY ##pochodne

    def langreneInterpolation(self, x):
        total = 0
        n = len(self.points)
        for i in range(n):
            xi, yi = self.points[i]
            Li = 1
            for j in range(n):
                if i != j:
                    xj, _ = self.points[j]
                    Li *= (x - xj) / (xi - xj)
            total += yi * Li
        return total

    def computeDifferenceQuotients(self):
        n = len(self.points)
        table = [[0] * n for _ in range(n)]
        for i in range(n):
            table[i][0] = self.points[i][1]

        for j in range(1, n):
            for i in range(n - j):
                numerator = table[i + 1][j - 1] - table[i][j - 1]
                denominator = self.points[i + j][0] - self.points[i][0]
                table[i][j] = numerator / denominator

        return [table[0][j] for j in range(n)]

    def newtonInterpolation(self, x):
        n = len(self.points)
        result = self.diffQuotients[0]
        productTerm = 1

        for i in range(1, n):
            productTerm *= x - self.points[i - 1][0]
            result += self.diffQuotients[i] * productTerm

        return result

    def computeProgressiveDifferences(self):
        n = len(self.points)
        diffs = [y for _, y in self.points]
        progressiveDiffs = [diffs]

        for i in range(1, n):
            currentDiffs = []
            for j in range(n - i):
                diff = progressiveDiffs[i - 1][j + 1] - progressiveDiffs[i - 1][j]
                currentDiffs.append(diff)
            progressiveDiffs.append(currentDiffs)

        return progressiveDiffs

    def newtonProgressiveInterpolation(self, x):
        n = len(self.points)
        result = self.progressiveDiffs[0][0]
        productTerm = 1

        for i in range(1, n):
            productTerm *= x - self.points[i - 1][0]
            result += (
                (self.progressiveDiffs[i][0] / (self.h**i))
                / self.factorial(i)
                * productTerm
            )

        return result

    @staticmethod
    def factorial(n):
        if n == 0 or n == 1:
            return 1
        else:
            return n * Interpolation.factorial(n - 1)

    def cubicSplineInterpolation(self, x):
        xPkt = [p[0] for p in self.points]
        yPkt = [p[1] for p in self.points]
        n = len(xPkt)
        A = [[0 for _ in range(n + 2)] for _ in range(n + 2)]
        b = yPkt + self.derY
        for i in range(n):
            xi = xPkt[i]
            A[i][:4] = [1, xi, xi**2, xi**3]
            for j in range(1, i):
                A[i][j + 3] = (xi - xPkt[j]) ** 3
        for i in range(2):
            xi = self.derX[i]
            A[n + i][1:4] = [1, 2 * xi, 3 * xi**2]
            if i == 1:
                for j in range(1, n - 1):
                    A[n + i][j + 3] = 3 * (xi - xPkt[j]) ** 2

        xWyniki = GaussElimination(A, b)

        hi = 0
        for i in range(1, n):
            if x < xPkt[i]:
                hi = i
                break
        wartoscWPunkcie = (
            xWyniki[0] + xWyniki[1] * x + xWyniki[2] * x**2 + xWyniki[3] * x**3
        )
        for i in range(1, hi):
            wartoscWPunkcie += xWyniki[i + 3] * (x - xPkt[i]) ** 3
        return wartoscWPunkcie
