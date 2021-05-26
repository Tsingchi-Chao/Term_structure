
import pandas as pd
import numpy as np
from collections import OrderedDict
import math
from bootstrapping.Bootstrapping import Bootstrapping
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

class LatticeTree:
    """LatticeTree
    This class will calibrate the lattice tree and get the short rate r and theta.
    You should initiallize K,sigma,delta t,Rmn and ytm.
    Here Rmn and ytm are the same with that in the Boostraping and here we use Rmn and ytm to initiallize the Boostrapping.
    """
    def __init__(self,
                 K:float,
                 sigma:float,
                 delta_t:float,
                 Rmn:dict,
                 ytm:dict):
        """
        :parameter
        :param K: The K in the Hull-White model formula, it stands for the mean-reversion speed.
        :param sigma: The sigma in the Hull-white model formula.
        :param delta_t: The delta t.
        """
        self.K=K
        self.sigma=sigma
        self.delta_t=delta_t
        self.n=math.floor(30 / delta_t)  #The longest period is 30 years. n is the num of step.
        self.J=math.floor(1/(2*K*delta_t))+1
        self.delta_r = sigma * math.sqrt(3 * delta_t)
        self.Rmn=Rmn
        self.ytm=ytm

    def get_k(self,
              j:int)->int:
        """
        k is a funtion of j, and it depends on the relationship between j and J.
        :parameter
        :param j:The step.
        :return:The k value.
        """
        J=self.J
        if j==-J:
            return -J+1
        elif abs(j)<J:
            return j
        elif j==J:
            return J-1
        else:
            print('There is something wrong,please check it.')

    def get_eta(self,
                j:int,
                n:int)->float:
        """
        This function is to calculating eta of j,n.
        parameter
        param j:The state.
        param n:The step.
        """
        eta_j_n=((1-self.K * self.delta_t)*j-LatticeTree.get_k(self,j))*self.delta_r
        return eta_j_n

    def get_p(self,
              j:int,
              n:int,
              x:int)->float:
        """
        This function is calculating probability by j,n,x and x may be k-1,k or k+1.
        parameter
        param j:The state.
        param n:The step.
        param x:It may be k-1,k or k+1.
        """
        k=LatticeTree.get_k(self,j)
        sigma=self.sigma
        delta_r=self.delta_r
        if x==k+1:
            p=0.5*(((sigma**2)*self.delta_t+(LatticeTree.get_eta(self,j,n))**2)/(delta_r**2)+LatticeTree.get_eta(self,j,n)/delta_r)
            return p
        elif x==k-1:
            p=0.5*(((sigma**2)*self.delta_t+(LatticeTree.get_eta(self,j,n))**2)/(delta_r**2)-LatticeTree.get_eta(self,j,n)/delta_r)
            return p
        elif x==k:
            p=1-((sigma**2)*self.delta_t+(LatticeTree.get_eta(self,j,n))**2)/(delta_r**2)
            return p
        else:
            p=0
            return p

    def get_r_theta_dict(self)->tuple:
        """
        It calculates short rate r and theta by loop from 1 to n-1.
        :parameter
        :param Rmn:The same with that in the Boostrapping.
        :param ytm:The same with that in the Boostraping.Because here we will use Boostrapping so we have to do initiallizaion.
        :return:r and theta.Both are dictionary whose first key is the step and the second key is the state.
        """
        boostrap=Bootstrapping(self.Rmn,self.ytm)
        Q = OrderedDict()  #The Arrow-Debreu Price.
        Q[0] = OrderedDict()
        Q[0][0] = 1
        r = OrderedDict()  #The short rate.
        r[0] = OrderedDict()
        r[0][0] = 0.00206
        theta = OrderedDict()  # The theta.
        n = self.n
        J=self.J
        delta_r=self.delta_r
        for j in range(1, n): #loop
            Q[j] = OrderedDict()
            r[j] = OrderedDict()
            for m in range(-1*min(j,J),min(j,J)+1):
                Q[j][m] = 0
                for i in range(-1*min(j-1,J), min(j-1,J) + 1):
                    Q[j][m] += Q[j-1][i]*LatticeTree.get_p(self,i,j-1,m) * math.exp(-r[j-1][i] * self.delta_t)
            for i in range(-1*min(j, J),min(j,J) + 1):
                r[j][i]=(1-self.K*self.delta_t)*r[j-1][0]+i*delta_r
            s = 0
            for i in range(-1*min(j,J),min(j,J)+1):
                s += Q[j][i] * math.exp(-1*r[j][i]*self.delta_t)
            theta[j-1] = (1/(self.K*(self.delta_t**2))) * math.log(s/boostrap.get_discount_rate_function(((j+1) *self.delta_t / 0.5)))
            for i in range(-1*min(j,J),min(j, J)+1):
                r[j][i]=r[j][i]+self.K*theta[j-1]*self.delta_t
        return r,theta


    def get_r_data(self)->pd.DataFrame:
        """
        This funtion is to get the dataframe of short rate r and theta.
        :return: It's a dataframe whose index is the term and columns are j,which means the state.
        """
        r=LatticeTree.get_r_theta_dict(self)[0]
        J=self.J
        data=pd.DataFrame()
        for i in range(1,self.n-1):
            for j in range(-J,J+1):
                try:
                    data.loc[i*self.delta_t,j] = r[i][j]
                except:
                    continue
        return data

    def get_theta_data(self)->pd.DataFrame:
        """
        This funtion is to get the dataframe of short rate r and theta.
        :return: It'a a dataframe whose index is the term and column is theta.
        """
        theta = LatticeTree.get_r_theta_dict(self)[1]
        data= pd.DataFrame()
        for i in range(self.n-1):
            data.loc[i * self.delta_t, 'theta'] = theta[i]
        return data

    def draw_r_picture(self):
        """
        This funtion is used to draw the pricture of the short rate tree.
        """
        data=LatticeTree.get_r_data(self)
        fig = plt.figure(figsize=[10, 6])
        ax = fig.add_subplot(1, 1, 1)
        for i in range(-self.J, self.J + 1):
            if i == 0:
                ax.plot(data[i], "*", c='b')
            else:
                ax.plot(data[i], "-")
        ax.set_xlabel('year')
        ax.set_ylabel('short rate')
        ax.set_title('The Tree')
        plt.show()

    def draw_theta_picture(self):
        """
        This funtion is used to draw the pricture of the theta.
        """
        data=LatticeTree.get_theta_data(self)
        fig = plt.figure(figsize=[10, 6])
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(data['theta'])
        ax.set_xlabel('year')
        ax.set_ylabel('theta')
        ax.set_title('The Theta')
        plt.show()

