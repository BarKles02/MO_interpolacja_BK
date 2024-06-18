from Interpolacja import Interpolation

if __name__ == "__main__":
    points2 = [(-4, 1127), (-2, 81), (0, 3), (2, 77), (4, 1023)]  ## to jest przykład
    xPochodne = [-4, 4]
    yPochodne = [-1093,1003]


    object = Interpolation(points2, xPochodne, yPochodne)
    # result = object.langreneInterpolation(3)
    # print(f'Interpolacja lagrangea: {result}')
    # print(object.computeDifferenceQuotients())
    # print(object.newtonInterpolation(3))
    
    # print(object.newtonProgressiveInterpolation(3))
    print(object.cubicSplineInterpolation(3))
