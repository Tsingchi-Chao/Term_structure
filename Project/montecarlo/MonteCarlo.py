
import pandas as pd
import numpy as np
import math
import time
from montecarlo.Parameter import Parameter
from bootstrapping.Bootstrapping import Bootstrapping
import warnings
warnings.filterwarnings("ignore")

class MonteCarlo:
    """Monte_Carlo
    A monte carlo to pricing the option of bond under Ho-Lee model.
    """
    def __init__(self,
                 ytm_discount_data:pd.DataFrame,
                 T0:float=0.5,
                 Tn:float=2,
                 delta_t:float=0.5,
                 coupon_rate:float=0.1,
                 strike_price:float=0.2):
        """
        :parameter
        :param T0:The maturity of the option.
        :param Tn:The matturity of the bond.
        :param delta_t:1/2 means there will be 2 coupon each year for the bond and 1/3 means 3 times.
        :param coupon_rate:The coupon rate of the bond.
        :param strike_price:The strike price of the option.
        """
        self.ytm_discount_data=ytm_discount_data
        self.T0=T0
        self.Tn=Tn
        self.delta_t=delta_t
        self.coupon_rate=coupon_rate
        self.strike_price=strike_price

    def get_PB(self,
               Ti:float,
               t:float=0)->float:
        """
        :parameter
        :param Ti:依次进行折现，这是指已经循环到的时间。
        :param t:指的是当前时间，如果是从开始计算就是0.
        :param P_B_inti:The inital P_B value.参考PPT。
        :returns
        :return 得到的是对应Ti下的P(T0,Ti)/B(T0).
        """
        P_B=self.ytm_discount_data.loc[Ti,'discount_rate']
        n=int((self.T0-t)/self.delta_t)
        para=Parameter(n,0.035,0.12)
        delta_W=para.delta_W()
        t=0
        for i in range(0,n):
            P_B=P_B*math.exp(-0.5*(para.Sigma(t,Ti)**2)*self.delta_t+para.Sigma(t,Ti)*delta_W[i])
            t+=self.delta_t
        return P_B

    def get_payoff(self)->float:
        """
        :return:求得每次模拟的payoff。
        """
        n=int((self.Tn-self.T0)/self.delta_t)
        EV=0 #指的是绝对payoff，即如果小于0就按照小于0表示。
        for Ti in np.arange(self.T0+self.delta_t,self.Tn+self.delta_t,self.delta_t):  #这里的公式参考PPT
            P_B=MonteCarlo.get_PB(self,Ti)
            EV=EV+P_B*self.coupon_rate*self.delta_t
            if Ti==self.Tn:
                EV=EV+P_B-self.strike_price*MonteCarlo.get_PB(self,self.T0)
        payoff=max(0,EV) #将EV与0比较，取最大的即为到期的payoff
        return payoff

    def run_monte_carlo(self,n:int):
        """
        :parameter
        :param n: The num of times you want to do montecarlo.
        :return: The average payoff which is also the price of the option.
        """
        data=pd.DataFrame()
        payoff_sum=0
        for j in range(n):
            payoff_sum+=MonteCarlo.get_payoff(self)
        average_payoff=payoff_sum/n
        return average_payoff