#!/usr/bin/env python
# coding: utf-8

# In[62]:


import pandas as pd
import numpy as np
import math
import time
import matplotlib.pyplot as plt


# 参考term structure PPT,Chpter 03-Part1

# In[76]:


class parameter:
    """parameter
    This class is used to calculate the sgima and delta_W.
    """
    
    def __init__(self,n:int,mu:float=0,std:float=math.sqrt(2)):
        """
        Args:
            n(float):The number of simulated values.这里生成随机数采用了一次性生成n个随机数，主要可以提升计算的速度。
            mu(float):The expected value of delta_W,here W is a standard brownian motion,and delta_W is the change.
            std(float):The standard deviation of the delta_W.
        """
        self.n=n
        self.mu=mu
        self.std=std
    
        
    def sigma(self,t:float,T:float)->float:
        """
        returns:
               The sigama of the P(t,T) under the Ho_Lee model.这是Ho_Lee model中最重要的变量之一，注意这是P(t,T)中的随机项，不是f(t,T)中的,
               这里的Sigama对T求导即得到f(t,T)中的随机项。
        """
        return 0.02*(T-t)
    
    def delta_W(self)->list:
        """
        Args:
             The number of the simulated values.
        Returns:
               The list of simulated value of the delta_W.
        """
        np.random.seed()  #这里没有用时间做种子，因为循环速度太快，导致了一个种子对应好几轮，进而出现伪随机数。
        s = np.random.normal(self.mu, self.std, self.n)  #get the simulated values
        return s


# In[92]:


class Monte_Carlo:
    """Monte_Carlo
    A monte carlo to pricing the option of bond under Ho-Lee model.
    """
    def __init__(self,T0:float=5,Tn:float=9,delta_t:float=0.5,coupon_rate:float=0.05,strike_price:float=1.01):
        """
        Args:
            T0(float):The maturity of the option.
            Tn(float):The matturity of the bond.
            delta_t(float):1/2 means there will be 2 coupon each year for the bond and 1/3 means 3 times.
            coupon_rate(float):The coupon rate of the bond.
            strike_price(float):The strike price of the option.
        """
        self.T0=T0
        self.Tn=Tn
        self.delta_t=delta_t
        self.coupon_rate=coupon_rate
        self.strike_price=strike_price
        
    def get_PB(self,Ti:float,t:float=0,P_B_init:float=0.05)->float:
        """
        Args:
            Ti(float):依次进行折现，这是指已经循环到的时间。
            t(float):指的是当前时间，如果是从开始计算就是0.
            P_B_inti(float):The inital P_B value.参考PPT。
        Returns:
            得到的是对应Ti下的P(T0,Ti)/B(T0).
        """
        P_B=P_B_init
        n=int((self.T0-t)/self.delta_t)
        para=parameter(n)
        delta_W=para.delta_W()
        t=0
        for i in range(0,n):
            P_B=P_B*math.exp(-0.5*(para.sigma(t,Ti)**2)*self.delta_t+para.sigma(t,Ti)*delta_W[i])
            t+=self.delta_t
        return P_B
    
    def get_payoff(self)->float:
        """
        Returns:
             求得每次模拟的payoff。
        """
        n=int((self.Tn-self.T0)/self.delta_t)
        EV=0 #指的是绝对payoff，即如果小于0就按照小于0表示。
        for Ti in np.arange(self.T0+self.delta_t,self.Tn+self.delta_t,self.delta_t):  #这里的公式参考PPT
            P_B=Monte_Carlo.get_PB(self,Ti)
            EV=EV+P_B*self.coupon_rate*self.delta_t
            if Ti==self.Tn:
                EV=EV+P_B-self.strike_price*Monte_Carlo.get_PB(self,self.T0)  
        payoff=max(0,EV) #将EV与0比较，取最大的即为到期的payoff
        return payoff


# In[95]:


def main():
    data=pd.DataFrame()
    n=2500 #最多模拟的次数
    for i in range(40,n,40):  #记录下不同模拟次数得到的结果
        payoff_sum=0
        for j in range(i):
            MTC=Monte_Carlo()
            payoff_sum+=MTC.get_payoff()
            payoff=payoff_sum/i
            data.loc[i,'payoff']=payoff
    V0=100*data[data.index>20]['payoff'].mean() #因为题目中是100元的面值，最后乘100即得到题目中的条件
    print('The price of the call option should be around',V0) 
    
    data['payoff']=100*data['payoff']
    fig = plt.figure(figsize=[10,6])
    ax= fig.add_subplot(1, 1, 1)
    ax.plot(data['payoff'],c='r',label='payoff')
    ax.legend()
    ax.axhline(V0,c='b')
    ax.set_xlabel('Times of simulation')
    ax.set_ylabel('payoff')
    ax.set_title('Pricing the bond option under one-factor Ho-Lee model b')


# In[ ]:


if __name__=="__main__":
    main()

