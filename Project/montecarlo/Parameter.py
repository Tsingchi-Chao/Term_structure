
import pandas as pd
import numpy as np
import math
import time
import warnings
warnings.filterwarnings("ignore")

class Parameter:
    """parameter
    This class is used to calculate the sgima and delta_W.
    """
    def __init__(self,n:int,sigma:float,K:float,mu:float=0,std:float=math.sqrt(0.5)):
        """
        :parameter
        :param n:The number of simulated values.这里生成随机数采用了一次性生成n个随机数，主要可以提升计算的速度。
        :param sigma: The sigma of the f(t,T),这里是f(t,T)中的随机项，不是P(t,T)中的随机项。
        :param K:The K of the Hull-white model.
        :param mu:The expected value of delta_W,here W is a standard brownian motion,and delta_W is the change.
        :param std:The standard deviation of the delta_W.
        """
        self.n=n
        self.K=K
        self.sigma=sigma
        self.mu=mu
        self.std=std

    def Sigma(self,t:float,T:float)->float:
        """
        returns:The sigama of the P(t,T) under the Hull-white model.这是Hull-white model中最重要的变量之一，注意这是P(t,T)中的随机项，不是f(t,T)中的,
                 这里的Sigama对T求导即得到f(t,T)中的随机项。
        """
        return self.sigma*(math.exp(-self.K*(T-t))-1)/self.K

    def delta_W(self)->list:
        """
        :parameter
        :param  The number of the simulated values.
        :return The list of simulated value of the delta_W.
        """
        np.random.seed()  #这里没有用时间做种子，因为循环速度太快，导致了一个种子对应好几轮，进而出现伪随机数。
        s = np.random.normal(self.mu, self.std, self.n)  #get the simulated values
        return s