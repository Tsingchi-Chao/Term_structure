
import pandas as pd
import numpy as np
from collections import OrderedDict
import math
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

class Bootstrapping:
    """Bootstrapping
    This class uses interest swap rate to bootstrap and get the zero-coupon yield curve and discount curve.
    It uses linear extrapolation.
    It uses bisection method to slove yield to maturity.
    In this case, we just need to confirm the first 0.5-year ytm and the interest swap rate dictionary then
    the class will calculate the whole ytm and the function of ytm and discount funtion.
    """
    def __init__(self,
                 Rmn:dict,
                 ytm:dict):
        """
        :parameter
        :param Rmn: It's a dictionary whose key is the term of the interest swap and the value is the interest swap rate.
                    And in this case we consider all the swaps are newly issued.
        :param ytm:It's a dictionary whose key is the the num of 0.5-year and value is the yield to maturity, we just need to
                   confirm the first 0.5-year yield to maturity.
        """
        self.Rmn=Rmn
        self.ytm=ytm

    def P_Tn(self,
             yn:float,
             n:int) -> float:
        """
        To calculate the present value.
        :parameter
        :param yn: The YTM for the zero-coupon bond.
        :param n: The number of the 0.5-year. It's very important for this.
        """
        return 1 / ((1 + yn * 0.5) **n)

    def f(self,
          yn: float,
          n: int,
          Rmn: dict,
          ytm: dict) -> float:
        """
        If we want to get yield to maturity we must make solve the equation and the analytical solution is a little
        complex so we decide to use bisection methed to solve the equation.
        And this is the funtion we want to make it equal 0 and if we solve it we can get Yield to maturity.
        :parameter
        :param yn: The ytm of the 0.5*n-year zero-coupon bond.
        :param n: The number of the 0.5-year.
        :param Rmn:The dictionary for the interest swap.
        :param ytm:The dictionary for the yield to maturity.
        """
        numerator = 1 - Bootstrapping.P_Tn(self,yn, n)
        denominator = 0
        for i in ytm.keys():
            denominator += 0.5 * Bootstrapping.P_Tn(self,ytm[i], i)
        ym = ytm[list(ytm.keys())[-1]]
        m = list(ytm.keys())[-1]
        for i in range(list(ytm.keys())[-1] + 1, n + 1):
            y = ym + (yn - ym) * (i - m) / (n - m)
            denominator += 0.5 * Bootstrapping.P_Tn(self,y, i)
        return numerator / denominator - self.Rmn[int(n / 2)]

    def bisection(self,
                  low_bound: float,
                  up_bound: float,
                  epsilon: float,
                  n: int,
                  Rmn: dict,
                  ytm: dict) -> float:
        """
        This is the bisection and by setting the low bound,up bound and the error we can solve the equaiton and get the yield to maturity.
        parameters
        param low_bound: The low bound you set.
        param up_bound: The up bound you set.
        param epsilon: The error.
        param n: The number of the 0.5-year.
        param Rmn:The dictionary for the interest swap.
        param ytm:The dictionary for the yield to maturity that have been solved.
        """
        for i in range(100):
            y = (low_bound + up_bound) / 2
            if abs(Bootstrapping.f(self,y, n, Rmn, ytm)) < epsilon:
                ym = ytm[list(ytm.keys())[-1]]
                m = list(ytm.keys())[-1]
                for j in range(list(ytm.keys())[-1] + 1, n + 1):
                    ytm[j] = ym + (y - ym) * (j - m) / (n - m)
                return ytm
            else:
                if Bootstrapping.f(self,y, n, Rmn, ytm) * Bootstrapping.f(self,up_bound, n, Rmn, ytm) < 0:
                    low_bound = y
                elif Bootstrapping.f(self,y, n, Rmn, ytm) * Bootstrapping.f(self,low_bound, n, Rmn, ytm) < 0:
                    up_bound = y

    def get_ytm_dict(self):
        """
        This function will get the whole dictionary of ytm. Here we set the low bound 0.001 and the high boung 0.1, the episilon 1e-10.
        :return:The dictionary of ytm whose key is 1,2,...,60 and value is the yield to maturity.
        """
        ytm=self.ytm
        for term in self.Rmn.keys():
            ytm = Bootstrapping.bisection(self,0.001, 0.1, 1e-10, 2 * term, self.Rmn, ytm)
        return ytm

    def get_ytm_discount_data(self):
        """
        This funtion is to get the dataframe of yield to maturity and discount rate value.
        :return: It'a dataframe whose index is the term and columns are yield to maturity and discount rate value.
        """
        ytm=Bootstrapping.get_ytm_dict(self)
        data = pd.DataFrame()
        for i in ytm.keys():
            data.loc[i / 2, 'Yield to maturity'] = ytm[i]
            data.loc[i / 2, 'discount_rate'] = Bootstrapping.P_Tn(self,ytm[i], i)
        return data

    def get_ytm_function(self,x: float):
        """
        Because get_ytm_dict will get a dictionary and get_ytm_discount_data
        will get dataframe but we want a function that we can input the term and then it give the yield to maturity.
        So this is that function.
        :parameter
        :param x:The number of 0.5-year.
        """
        ytm = Bootstrapping.get_ytm_dict(self)
        if x<1:
            return ytm[1]
        if x>=1 and x<=4:
            return ytm[1]+((ytm[4]-ytm[1])/(4-1))*(x-1)
        elif x>4 and x<=6:
            return ytm[4]+((ytm[6]-ytm[4])/(6-4))*(x-4)
        elif x>6 and x<=8:
            return ytm[6]+((ytm[8]-ytm[6])/(8-6))*(x-6)
        elif x>8 and x<=10:
            return ytm[8]+((ytm[10]-ytm[8])/(10-8))*(x-8)
        elif x>10 and x<=14:
            return ytm[10]+((ytm[14]-ytm[10])/(14-10))*(x-10)
        elif x>14 and x<=20:
            return ytm[14]+((ytm[20]-ytm[14])/(20-14))*(x-14)
        elif x>20 and x<=40:
            return ytm[20]+((ytm[40]-ytm[20])/(40-20))*(x-20)
        elif x>40 and x<=60:
            return ytm[40]+((ytm[60]-ytm[40])/(60-40))*(x-40)

    def get_discount_rate_function(self,x: float):
        """
        Because get_ytm_discount_data will get dataframe but we want a function that we can
        input the term and then it give the discount rate.
        So this is that function.
        parameter
        param x:The number of 0.5-year.
        """
        discount_rate = Bootstrapping.P_Tn(self,Bootstrapping.get_ytm_function(self,x), x)
        return discount_rate

    def draw_yield_curve(self):
        """
        This function is used to draw the zero-coupon bond yield curve.
        """
        data = Bootstrapping.get_ytm_discount_data(self)
        fig = plt.figure(figsize=[10, 6])
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(data['Yield to maturity'])
        ax.set_xlabel('year')
        ax.set_ylabel('rate')
        ax.set_title('Zero-coupon yield curve')
        plt.show()

    def draw_discount_curve(self):
        """
        This function is used to draw the zero-coupon bond yield curve.
        """
        data=Bootstrapping.get_ytm_discount_data(self)
        fig = plt.figure(figsize=[10, 6])
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(data['discount_rate'])
        ax.set_xlabel('Term')
        ax.set_ylabel('value')
        ax.set_title('Discount Curves')
        plt.show()



